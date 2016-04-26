__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import inspect
import os
import traceback
import collections
from urllib.request import urlopen
from urllib.error import URLError

from snakemake.utils import format
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError

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


def script(basedir, path, input, output, params, wildcards, threads, resources,
           log, config):
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
                try:
                    exec(compile(source.read().decode(), path, "exec"), {
                        "snakemake": Snakemake(input, output, params, wildcards,
                                               threads, resources, log, config)
                    })
                except (Exception, BaseException) as ex:
                    raise WorkflowError("".join(traceback.format_exception(type(ex), ex, ex.__traceback__)))
            elif path.endswith(".R"):
                try:
                    import rpy2.robjects as robjects
                except ImportError:
                    raise ValueError(
                        "Python 3 package rpy2 needs to be installed to use the R function.")
                with urlopen(path) as source:
                    preamble = """
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
                    """.format(REncoder.encode_namedlist(input),
                               REncoder.encode_namedlist(output),
                               REncoder.encode_namedlist(params),
                               REncoder.encode_namedlist(wildcards), threads,
                               REncoder.encode_namedlist(log),
                               REncoder.encode_namedlist({
                                   name: value
                                   for name, value in resources.items()
                                   if name != "_cores" and name != "_nodes"
                               }), REncoder.encode_dict(config))
                    logger.debug(preamble)
                    source = preamble + source.read().decode()
                    robjects.r(source)
            else:
                raise ValueError(
                    "Unsupported script: Expecting either Python (.py) or R (.R) script.")
    except URLError as e:
        raise WorkflowError(e)
