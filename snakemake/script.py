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
import subprocess
import collections
import re
from urllib.request import urlopen, pathname2url
from urllib.error import URLError
from itertools import islice

from snakemake.utils import format
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.common import MIN_PY_VERSION, escape_backslash, SNAKEMAKE_SEARCHPATH
from snakemake.io import git_content, split_git_path
from snakemake import singularity


PY_VER_RE = re.compile("Python (?P<ver_min>\d+\.\d+).*")
# TODO use this to find the right place for inserting the preamble
PY_PREAMBLE_RE = re.compile(r"from( )+__future__( )+import.*?(?P<end>[;\n])")


class REncoder:
    """Encoding Pyton data structures into R."""

    @classmethod
    def encode_numeric(cls, value):
        if value is None:
            return "as.numeric(NA)"
        return str(value)

    @classmethod
    def encode_value(cls, value):
        if value is None:
            return "NULL"
        elif isinstance(value, str):
            return repr(value)
        elif isinstance(value, dict):
            return cls.encode_dict(value)
        elif isinstance(value, bool):
            return "TRUE" if value else "FALSE"
        elif isinstance(value, int) or isinstance(value, float):
            return str(value)
        elif isinstance(value, collections.Iterable):
            # convert all iterables to vectors
            return cls.encode_list(value)
        else:
            # Try to convert from numpy if numpy is present
            try:
                import numpy as np

                if isinstance(value, np.number):
                    return str(value)
            except ImportError:
                pass
        raise ValueError("Unsupported value for conversion into R: {}".format(value))

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
        positional = ", ".join(map(cls.encode_value, namedlist))
        named = cls.encode_items(namedlist.items())
        source = "list("
        if positional:
            source += positional
        if named:
            source += ", " + named
        source += ")"
        return source


class JuliaEncoder:
    """Encoding Pyton data structures into Julia."""

    @classmethod
    def encode_value(cls, value):
        if value is None:
            return "nothing"
        elif isinstance(value, str):
            return repr(value)
        elif isinstance(value, dict):
            return cls.encode_dict(value)
        elif isinstance(value, bool):
            return "true" if value else "false"
        elif isinstance(value, int) or isinstance(value, float):
            return str(value)
        elif isinstance(value, collections.Iterable):
            # convert all iterables to vectors
            return cls.encode_list(value)
        else:
            # Try to convert from numpy if numpy is present
            try:
                import numpy as np

                if isinstance(value, np.number):
                    return str(value)
            except ImportError:
                pass
        raise ValueError(
            "Unsupported value for conversion into Julia: {}".format(value)
        )

    @classmethod
    def encode_list(cls, l):
        return "[{}]".format(", ".join(map(cls.encode_value, l)))

    @classmethod
    def encode_items(cls, items):
        def encode_item(item):
            name, value = item
            return '"{}" => {}'.format(name, cls.encode_value(value))

        return ", ".join(map(encode_item, items))

    @classmethod
    def encode_positional_items(cls, namedlist):
        encoded = ""
        for index, value in enumerate(namedlist):
            encoded += "{} => {}, ".format(index + 1, cls.encode_value(value))
        return encoded

    @classmethod
    def encode_dict(cls, d):
        d = "Dict({})".format(cls.encode_items(d.items()))
        return d

    @classmethod
    def encode_namedlist(cls, namedlist):
        positional = cls.encode_positional_items(namedlist)
        named = cls.encode_items(namedlist.items())
        source = "Dict("
        if positional:
            source += positional
        if named:
            source += named
        source += ")"
        return source


class Snakemake:
    def __init__(
        self,
        input,
        output,
        params,
        wildcards,
        threads,
        resources,
        log,
        config,
        rulename,
        bench_iteration,
        scriptdir=None,
    ):
        # convert input and output to plain strings as some remote objects cannot
        # be pickled
        self.input = input.plainstrings()
        self.output = output.plainstrings()
        self.params = params
        self.wildcards = wildcards
        self.threads = threads
        self.resources = resources
        self.log = log.plainstrings()
        self.config = config
        self.rule = rulename
        self.bench_iteration = bench_iteration
        self.scriptdir = scriptdir

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


def get_source(path, basedir="."):
    source = None
    if not path.startswith("http") and not path.startswith("git+file"):
        if path.startswith("file://"):
            path = path[7:]
        elif path.startswith("file:"):
            path = path[5:]
        if not os.path.isabs(path):
            path = os.path.abspath(os.path.join(basedir, path))
        path = "file://" + path
    path = format(path, stepout=1)
    if path.startswith("file://"):
        sourceurl = "file:" + pathname2url(path[7:])
    elif path.startswith("git+file"):
        source = git_content(path)
        (root_path, file_path, version) = split_git_path(path)
        path = path.rstrip("@" + version)
    else:
        sourceurl = path

    language = None
    if path.endswith(".py"):
        language = "python"
    elif path.endswith(".R"):
        language = "r"
    elif path.endswith(".Rmd"):
        language = "rmarkdown"
    elif path.endswith(".jl"):
        language = "julia"

    if source is None:
        with urlopen(sourceurl) as source:
            return path, source.read(), language
    else:
        return path, source, language


def script(
    path,
    basedir,
    input,
    output,
    params,
    wildcards,
    threads,
    resources,
    log,
    config,
    rulename,
    conda_env,
    singularity_img,
    singularity_args,
    bench_record,
    jobid,
    bench_iteration,
    cleanup_scripts,
    shadow_dir,
):
    """
    Load a script from the given basedir + path and execute it.
    Supports Python 3 and R.
    """

    f = None
    try:
        path, source, language = get_source(path, basedir)
        if language == "python":
            wrapper_path = path[7:] if path.startswith("file://") else path
            snakemake = Snakemake(
                input,
                output,
                params,
                wildcards,
                threads,
                resources,
                log,
                config,
                rulename,
                bench_iteration,
                os.path.dirname(wrapper_path),
            )
            snakemake = pickle.dumps(snakemake)
            # Obtain search path for current snakemake module.
            # The module is needed for unpickling in the script.
            # We append it at the end (as a fallback).
            searchpath = SNAKEMAKE_SEARCHPATH
            if singularity_img is not None:
                searchpath = singularity.SNAKEMAKE_MOUNTPOINT
            searchpath = '"{}"'.format(searchpath)
            # For local scripts, add their location to the path in case they use path-based imports
            if path.startswith("file://"):
                searchpath += ', "{}"'.format(os.path.dirname(path[7:]))
            preamble = textwrap.dedent(
                """
            ######## Snakemake header ########
            import sys; sys.path.extend([{searchpath}]); import pickle; snakemake = pickle.loads({snakemake}); from snakemake.logging import logger; logger.printshellcmds = {printshellcmds}; __real_file__ = __file__; __file__ = {file_override};
            ######## Original script #########
            """
            ).format(
                searchpath=escape_backslash(searchpath),
                snakemake=snakemake,
                printshellcmds=logger.printshellcmds,
                file_override=repr(os.path.realpath(wrapper_path)),
            )
        elif language == "r" or language == "rmarkdown":
            preamble = textwrap.dedent(
                """
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
                    config = "list",
                    rule = "character",
                    bench_iteration = "numeric",
                    scriptdir = "character",
                    source = "function"
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
                config = {},
                rule = {},
                bench_iteration = {},
                scriptdir = {},
                source = function(...){{
                    wd <- getwd()
                    setwd(snakemake@scriptdir)
                    source(...)
                    setwd(wd)
                }}
            )

            ######## Original script #########
            """
            ).format(
                REncoder.encode_namedlist(input),
                REncoder.encode_namedlist(output),
                REncoder.encode_namedlist(params),
                REncoder.encode_namedlist(wildcards),
                threads,
                REncoder.encode_namedlist(log),
                REncoder.encode_namedlist(
                    {
                        name: value
                        for name, value in resources.items()
                        if name != "_cores" and name != "_nodes"
                    }
                ),
                REncoder.encode_dict(config),
                REncoder.encode_value(rulename),
                REncoder.encode_numeric(bench_iteration),
                REncoder.encode_value(
                    os.path.dirname(path[7:])
                    if path.startswith("file://")
                    else os.path.dirname(path)
                ),
            )
        elif language == "julia":
            preamble = textwrap.dedent(
                """
                    ######## Snakemake header ########
                    struct Snakemake
                        input::Dict
                        output::Dict
                        params::Dict
                        wildcards::Dict
                        threads::Int64
                        log::Dict
                        resources::Dict
                        config::Dict
                        rule::String
                        bench_iteration
                        scriptdir::String
                        #source::Any
                    end
                    snakemake = Snakemake(
                        {}, #input::Dict
                        {}, #output::Dict
                        {}, #params::Dict
                        {}, #wildcards::Dict
                        {}, #threads::Int64
                        {}, #log::Dict
                        {}, #resources::Dict
                        {}, #config::Dict
                        {}, #rule::String
                        {}, #bench_iteration::Int64
                        {}, #scriptdir::String
                        #, #source::Any
                    )
                    ######## Original script #########
                    """.format(
                    JuliaEncoder.encode_namedlist(input),
                    JuliaEncoder.encode_namedlist(output),
                    JuliaEncoder.encode_namedlist(params),
                    JuliaEncoder.encode_namedlist(wildcards),
                    JuliaEncoder.encode_value(threads),
                    JuliaEncoder.encode_namedlist(log),
                    JuliaEncoder.encode_namedlist(
                        {
                            name: value
                            for name, value in resources.items()
                            if name != "_cores" and name != "_nodes"
                        }
                    ),
                    JuliaEncoder.encode_dict(config),
                    JuliaEncoder.encode_value(rulename),
                    JuliaEncoder.encode_value(bench_iteration),
                    JuliaEncoder.encode_value(
                        os.path.dirname(path[7:])
                        if path.startswith("file://")
                        else os.path.dirname(path)
                    ),
                ).replace(
                    "'", '"'
                )
            )
        else:
            raise ValueError(
                "Unsupported script: Expecting either Python (.py), R (.R), RMarkdown (.Rmd) or Julia (.jl) script."
            )

        dir = ".snakemake/scripts"
        os.makedirs(dir, exist_ok=True)

        with tempfile.NamedTemporaryFile(
            suffix="." + os.path.basename(path), dir=dir, delete=False
        ) as f:
            if not language == "rmarkdown":
                f.write(preamble.encode())
                f.write(source)
            else:
                # Insert Snakemake object after the RMarkdown header
                code = source.decode()
                pos = next(islice(re.finditer(r"---\n", code), 1, 2)).start() + 3
                f.write(str.encode(code[:pos]))
                preamble = textwrap.dedent(
                    """
                    ```{r, echo=FALSE, message=FALSE, warning=FALSE}
                    %s
                    ```
                    """
                    % preamble
                )
                f.write(preamble.encode())
                f.write(str.encode(code[pos:]))

        if language == "python":
            py_exec = sys.executable
            if conda_env is not None:
                py = os.path.join(conda_env, "bin", "python")
                if os.path.exists(py):
                    out = subprocess.check_output(
                        [py, "--version"],
                        stderr=subprocess.STDOUT,
                        universal_newlines=True,
                    )
                    ver = tuple(
                        map(int, PY_VER_RE.match(out).group("ver_min").split("."))
                    )
                    if ver >= MIN_PY_VERSION:
                        # Python version is new enough, make use of environment
                        # to execute script
                        py_exec = "python"
                    else:
                        logger.warning(
                            "Conda environment defines Python "
                            "version < {0}.{1}. Using Python of the "
                            "master process to execute "
                            "script. Note that this cannot be avoided, "
                            "because the script uses data structures from "
                            "Snakemake which are Python >={0}.{1} "
                            "only.".format(*MIN_PY_VERSION)
                        )
            if singularity_img is not None:
                # use python from image
                py_exec = "python"
            # use the same Python as the running process or the one from the environment
            shell("{py_exec} {f.name:q}", bench_record=bench_record)
        elif language == "r":
            if conda_env is not None and "R_LIBS" in os.environ:
                logger.warning(
                    "R script job uses conda environment but "
                    "R_LIBS environment variable is set. This "
                    "is likely not intended, as R_LIBS can "
                    "interfere with R packages deployed via "
                    "conda. Consider running `unset R_LIBS` or "
                    "remove it entirely before executing "
                    "Snakemake."
                )
            shell("Rscript --vanilla {f.name:q}", bench_record=bench_record)
        elif language == "rmarkdown":
            if len(output) != 1:
                raise WorkflowError(
                    "RMarkdown scripts (.Rmd) may only have a single output file."
                )
            out = os.path.abspath(output[0])
            shell(
                'Rscript --vanilla -e \'rmarkdown::render("{f.name}", output_file="{out}", quiet=TRUE, knit_root_dir = "{workdir}", params = list(rmd="{f.name}"))\'',
                bench_record=bench_record,
                workdir=os.getcwd(),
            )
        elif language == "julia":
            shell("julia {f.name:q}", bench_record=bench_record)

    except URLError as e:
        raise WorkflowError(e)
    finally:
        if f and cleanup_scripts:
            os.remove(f.name)
        else:
            logger.debug("Not cleaning up %s" % f.name)
