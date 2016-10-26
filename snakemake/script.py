__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import inspect
import os
import tempfile
import textwrap
import sys
import pickle
import traceback
import collections
from urllib.request import urlopen
from urllib.error import URLError

from snakemake.utils import format
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
from snakemake.shell import shell

class REncoder:
    """Encoding Pyton data structures into R."""

    @classmethod
    def encode_value(cls, value):
        if isinstance(value, str):
            return repr(value)
        if isinstance(value, collections.Iterable):
            # convert all iterables to vectors
            return cls.encode_list(value)
        elif isinstance(value, dict):
            return cls.encode_dict(value)
        elif isinstance(value, bool):
            return "TRUE" if value else "FALSE"
        elif isinstance(value, int) or isinstance(value, float):
            return str(value)
        else:
            # Try to convert from numpy if numpy is present
            try:
                import numpy as np
                if isinstance(value, np.number):
                    return str(value)
            except ImportError:
                pass
        raise ValueError(
            "Unsupported value for conversion into R: {}".format(value))

    @classmethod
    def encode_list(cls, l):
        return "c({})".format(", ".join(map(cls.encode_value, l)))

    @classmethod
    def encode_items(cls, items):
        def encode_item(item):
            name, value = item
            return '"{}" = {}'.format(name, cls.encode_value(value))

        return ", ".join(map(encode_item, items))

    @classmethod
    def encode_dict(cls, d):
        d = "list({})".format(cls.encode_items(d.items()))
        return d

    @classmethod
    def encode_namedlist(cls, namedlist):
        positional = cls.encode_list(namedlist)
        named = cls.encode_items(namedlist.items())
        source = "list("
        if positional != "c()":
            source += positional
        if named:
            source += ", " + named
        source += ")"
        return source


class Snakemake:
    def __init__(self, input, output, params, wildcards, threads, resources,
                 log, config):
        self.input = input
        self.output = output
        self.params = params
        self.wildcards = wildcards
        self.threads = threads
        self.resources = resources
        self.log = log
        self.config = config

    def log_fmt_shell(self, stdout=True, stderr=True, append=False):
        """
        Return a shell redirection string to be used in `shell()` calls

        This function allows scripts and wrappers support optional `log` files
        specified in the calling rule.  If no `log` was specified, then an
        empty string "" is returned, regardless of the values of `stdout`,
        `stderr`, and `append`.

        Parameters
        ---------

        stdout : bool
            Send stdout to log

        stderr : bool
            Send stderr to log

        append : bool
            Do not overwrite the log file. Useful for sending output of
            multiple commands to the same log. Note however that the log will
            not be truncated at the start.

        The following table describes the output:

        -------- -------- -------- ----- -------------
        stdout   stderr   append   log   return value
        -------- -------- -------- ----- ------------
        True     True     True     fn    >> fn 2>&1
        True     False    True     fn    >> fn
        False    True     True     fn    2>> fn
        True     True     False    fn    > fn 2>&1
        True     False    False    fn    > fn
        False    True     False    fn    2> fn
        any      any      any      None  ""
        -------- -------- -------- ----- -----------
        """
        if not self.log:
            return ""
        lookup = {
            (True, True, True): " >> {0} 2>&1",
            (True, False, True): " >> {0}",
            (False, True, True): " 2>> {0}",
            (True, True, False): " > {0} 2>&1",
            (True, False, False): " > {0}",
            (False, True, False): " 2> {0}",
        }
        return lookup[(stdout, stderr, append)].format(self.log)


def script(path, basedir, input, output, params, wildcards, threads, resources,
           log, config, conda_env):
    """
    Load a script from the given basedir + path and execute it.
    Supports Python 3 and R.
    """
    if not path.startswith("http"):
        if path.startswith("file://"):
            path = path[7:]
        if not os.path.isabs(path):
            path = os.path.abspath(os.path.join(basedir, path))
        path = "file://" + path
    path = format(path, stepout=1)

    try:
        with urlopen(path) as source:
            if path.endswith(".py"):
                snakemake = Snakemake(input, output, params, wildcards,
                                      threads, resources, log, config)
                snakemake = pickle.dumps(snakemake)
                # obtain search path for current snakemake module
                # the module is needed for unpickling in the script
                snakemake_path = os.path.dirname(os.path.dirname(__file__))
                preamble = textwrap.dedent("""
                ######## Snakemake header ########
                import sys; sys.path.insert(0, "{}"); sys.path.extend({}); import pickle; snakemake = pickle.loads({})
                ######## Original script #########
                """).format(snakemake_path, sys.path, snakemake)
            elif path.endswith(".R"):
                preamble = textwrap.dedent("""
                ######## Snakemake header ########
                library(methods)
                Snakemake <- setClass(
                    "Snakemake",
                    slots = c(
                        input = "list",
                        output = "list",
                        params = "list",
                        wildcards = "list",
                        threads = "numeric",
                        log = "list",
                        resources = "list",
                        config = "list"
                    )
                )
                snakemake <- Snakemake(
                    input = {},
                    output = {},
                    params = {},
                    wildcards = {},
                    threads = {},
                    log = {},
                    resources = {},
                    config = {}
                )
                ######## Original script #########
                """).format(REncoder.encode_namedlist(input),
                           REncoder.encode_namedlist(output),
                           REncoder.encode_namedlist(params),
                           REncoder.encode_namedlist(wildcards), threads,
                           REncoder.encode_namedlist(log),
                           REncoder.encode_namedlist({
                               name: value
                               for name, value in resources.items()
                               if name != "_cores" and name != "_nodes"
                           }), REncoder.encode_dict(config))
            else:
                raise ValueError(
                    "Unsupported script: Expecting either Python (.py) or R (.R) script.")

            dir = ".snakemake/scripts"
            os.makedirs(dir, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                suffix="." + os.path.basename(path),
                prefix="",
                dir=dir,
                delete=False) as f:
                f.write(preamble.encode())
                f.write(source.read())
            if path.endswith(".py"):
                shell("{python} {f.name}",
                      python=sys.executable) # always use the same Python as the running process
            elif path.endswith(".R"):
                shell("Rscript {f.name}")
            os.remove(f.name)

    except URLError as e:
        raise WorkflowError(e)
