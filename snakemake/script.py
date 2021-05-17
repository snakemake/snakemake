__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import inspect
import itertools
import os
import tempfile
import textwrap
import sys
import pickle
import subprocess
import collections
import re
from abc import ABC, abstractmethod
from urllib.request import urlopen, pathname2url
from urllib.error import URLError

from snakemake.utils import format
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.common import (
    MIN_PY_VERSION,
    SNAKEMAKE_SEARCHPATH,
    ON_WINDOWS,
    smart_join,
    is_local_file,
)
from snakemake.io import git_content, split_git_path
from snakemake.deployment import singularity


# TODO use this to find the right place for inserting the preamble
PY_PREAMBLE_RE = re.compile(r"from( )+__future__( )+import.*?(?P<end>[;\n])")


class Snakemake:
    def __init__(
        self,
        input_,
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
        self.input = input_._plainstrings()
        self.output = output._plainstrings()
        self.params = params
        self.wildcards = wildcards
        self.threads = threads
        self.resources = resources
        self.log = log._plainstrings()
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
        elif isinstance(value, collections.abc.Iterable):
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
        elif isinstance(value, collections.abc.Iterable):
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


class ScriptBase(ABC):
    editable = False

    def __init__(
        self,
        path,
        source,
        basedir,
        input_,
        output,
        params,
        wildcards,
        threads,
        resources,
        log,
        config,
        rulename,
        conda_env,
        container_img,
        singularity_args,
        env_modules,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
    ):
        self.path = path
        self.source = source

        self.basedir = basedir
        self.input = input_
        self.output = output
        self.params = params
        self.wildcards = wildcards
        self.threads = threads
        self.resources = resources
        self.log = log
        self.config = config
        self.rulename = rulename
        self.conda_env = conda_env
        self.container_img = container_img
        self.singularity_args = singularity_args
        self.env_modules = env_modules
        self.bench_record = bench_record
        self.jobid = jobid
        self.bench_iteration = bench_iteration
        self.cleanup_scripts = cleanup_scripts
        self.shadow_dir = shadow_dir

    def evaluate(self, edit=False):
        assert not edit or self.editable

        fd = None
        try:
            # generate preamble
            preamble = self.get_preamble()

            # write script
            dir_ = ".snakemake/scripts"
            os.makedirs(dir_, exist_ok=True)

            with tempfile.NamedTemporaryFile(
                suffix="." + os.path.basename(self.path), dir=dir_, delete=False
            ) as fd:
                self.write_script(preamble, fd)

            # execute script
            self.execute_script(fd.name, edit=edit)
        except URLError as e:
            raise WorkflowError(e)
        finally:
            if fd and self.cleanup_scripts:
                os.remove(fd.name)
            else:
                if fd:
                    logger.warning("Not cleaning up %s" % fd.name)
                else:
                    # nothing to clean up (TODO: ??)
                    pass

    @property
    def local_path(self):
        path = self.path[7:]
        if not os.path.isabs(path):
            return smart_join(self.basedir, path)
        return path

    @abstractmethod
    def get_preamble(self):
        ...

    @abstractmethod
    def write_script(self, preamble, fd):
        ...

    @abstractmethod
    def execute_script(self, fname, edit=False):
        ...

    def _execute_cmd(self, cmd, **kwargs):
        return shell(
            cmd,
            bench_record=self.bench_record,
            conda_env=self.conda_env,
            container_img=self.container_img,
            shadow_dir=self.shadow_dir,
            env_modules=self.env_modules,
            singularity_args=self.singularity_args,
            **kwargs
        )


class PythonScript(ScriptBase):
    @staticmethod
    def generate_preamble(
        path,
        source,
        basedir,
        input_,
        output,
        params,
        wildcards,
        threads,
        resources,
        log,
        config,
        rulename,
        conda_env,
        container_img,
        singularity_args,
        env_modules,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
        preamble_addendum="",
    ):
        wrapper_path = path[7:] if path.startswith("file://") else path
        snakemake = Snakemake(
            input_,
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
        if container_img is not None:
            searchpath = singularity.SNAKEMAKE_MOUNTPOINT
        searchpath = repr(searchpath)
        # For local scripts, add their location to the path in case they use path-based imports
        if path.startswith("file://"):
            searchpath += ", " + repr(os.path.dirname(path[7:]))

        return textwrap.dedent(
            """
        ######## snakemake preamble start (automatically inserted, do not edit) ########
        import sys; sys.path.extend([{searchpath}]); import pickle; snakemake = pickle.loads({snakemake}); from snakemake.logging import logger; logger.printshellcmds = {printshellcmds}; {preamble_addendum}
        ######## snakemake preamble end #########
        """
        ).format(
            searchpath=searchpath,
            snakemake=snakemake,
            printshellcmds=logger.printshellcmds,
            preamble_addendum=preamble_addendum,
        )

    def get_preamble(self):
        wrapper_path = self.path[7:] if self.path.startswith("file://") else self.path
        preamble_addendum = (
            "__real_file__ = __file__; __file__ = {file_override};".format(
                file_override=repr(os.path.realpath(wrapper_path))
            )
        )

        return PythonScript.generate_preamble(
            self.path,
            self.source,
            self.basedir,
            self.input,
            self.output,
            self.params,
            self.wildcards,
            self.threads,
            self.resources,
            self.log,
            self.config,
            self.rulename,
            self.conda_env,
            self.container_img,
            self.singularity_args,
            self.env_modules,
            self.bench_record,
            self.jobid,
            self.bench_iteration,
            self.cleanup_scripts,
            self.shadow_dir,
            preamble_addendum=preamble_addendum,
        )

    def write_script(self, preamble, fd):
        fd.write(preamble.encode())
        fd.write(self.source)

    def _is_python_env(self):
        if self.conda_env is not None:
            prefix = os.path.join(self.conda_env, "bin")
        elif self.env_modules is not None:
            prefix = self._execute_cmd("echo $PATH", read=True).split(":")[0]
        else:
            raise NotImplementedError()
        return os.path.exists(os.path.join(prefix, "python"))

    def _get_python_version(self):
        out = self._execute_cmd(
            "python -c \"import sys; print('.'.join(map(str, sys.version_info[:2])))\"",
            read=True,
        )
        return tuple(map(int, out.strip().split(".")))

    def execute_script(self, fname, edit=False):
        py_exec = sys.executable
        if self.container_img is not None:
            # use python from image
            py_exec = "python"
        elif self.conda_env is not None or self.env_modules is not None:
            if self._is_python_env():
                py_version = self._get_python_version()
                # If version is None, all fine, because host python usage is intended.
                if py_version is not None:
                    if py_version >= MIN_PY_VERSION:
                        # Python version is new enough, make use of environment
                        # to execute script
                        py_exec = "python"
                    else:
                        logger.warning(
                            "Environment defines Python "
                            "version < {0}.{1}. Using Python of the "
                            "main process to execute "
                            "script. Note that this cannot be avoided, "
                            "because the script uses data structures from "
                            "Snakemake which are Python >={0}.{1} "
                            "only.".format(*MIN_PY_VERSION)
                        )

        if ON_WINDOWS:
            # use forward slashes so script command still works even if
            # bash is configured as executable on Windows
            py_exec = py_exec.replace("\\", "/")
        # use the same Python as the running process or the one from the environment
        self._execute_cmd("{py_exec} {fname:q}", py_exec=py_exec, fname=fname)


class RScript(ScriptBase):
    @staticmethod
    def generate_preamble(
        path,
        source,
        basedir,
        input_,
        output,
        params,
        wildcards,
        threads,
        resources,
        log,
        config,
        rulename,
        conda_env,
        container_img,
        singularity_args,
        env_modules,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
        preamble_addendum="",
    ):
        return textwrap.dedent(
            """
        ######## snakemake preamble start (automatically inserted, do not edit) ########
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
        {preamble_addendum}

        ######## snakemake preamble end #########
        """
        ).format(
            REncoder.encode_namedlist(input_),
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
            preamble_addendum=preamble_addendum,
        )

    def get_preamble(self):
        return RScript.generate_preamble(
            self.path,
            self.source,
            self.basedir,
            self.input,
            self.output,
            self.params,
            self.wildcards,
            self.threads,
            self.resources,
            self.log,
            self.config,
            self.rulename,
            self.conda_env,
            self.container_img,
            self.singularity_args,
            self.env_modules,
            self.bench_record,
            self.jobid,
            self.bench_iteration,
            self.cleanup_scripts,
            self.shadow_dir,
        )

    def write_script(self, preamble, fd):
        fd.write(preamble.encode())
        fd.write(self.source)

    def execute_script(self, fname, edit=False):
        if self.conda_env is not None and "R_LIBS" in os.environ:
            logger.warning(
                "R script job uses conda environment but "
                "R_LIBS environment variable is set. This "
                "is likely not intended, as R_LIBS can "
                "interfere with R packages deployed via "
                "conda. Consider running `unset R_LIBS` or "
                "remove it entirely before executing "
                "Snakemake."
            )
        self._execute_cmd("Rscript --vanilla {fname:q}", fname=fname)


class RMarkdown(ScriptBase):
    def get_preamble(self):
        return textwrap.dedent(
            """
        ######## snakemake preamble start (automatically inserted, do not edit) ########
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

        ######## snakemake preamble end #########
        """
        ).format(
            REncoder.encode_namedlist(self.input),
            REncoder.encode_namedlist(self.output),
            REncoder.encode_namedlist(self.params),
            REncoder.encode_namedlist(self.wildcards),
            self.threads,
            REncoder.encode_namedlist(self.log),
            REncoder.encode_namedlist(
                {
                    name: value
                    for name, value in self.resources.items()
                    if name != "_cores" and name != "_nodes"
                }
            ),
            REncoder.encode_dict(self.config),
            REncoder.encode_value(self.rulename),
            REncoder.encode_numeric(self.bench_iteration),
            REncoder.encode_value(
                os.path.dirname(self.path[7:])
                if self.path.startswith("file://")
                else os.path.dirname(self.path)
            ),
        )

    def write_script(self, preamble, fd):
        # Insert Snakemake object after the RMarkdown header
        code = self.source.decode()
        pos = next(itertools.islice(re.finditer(r"---\n", code), 1, 2)).start() + 3
        fd.write(str.encode(code[:pos]))
        preamble = textwrap.dedent(
            """
            ```{r, echo=FALSE, message=FALSE, warning=FALSE}
            %s
            ```
            """
            % preamble
        )
        fd.write(preamble.encode())
        fd.write(str.encode(code[pos:]))

    def execute_script(self, fname, edit=False):
        if len(self.output) != 1:
            raise WorkflowError(
                "RMarkdown scripts (.Rmd) may only have a single output file."
            )
        out = os.path.abspath(self.output[0])
        self._execute_cmd(
            'Rscript --vanilla -e \'rmarkdown::render("{fname}", output_file="{out}", quiet=TRUE, knit_root_dir = "{workdir}", params = list(rmd="{fname}"))\'',
            fname=fname,
            out=out,
            workdir=os.getcwd(),
        )


class JuliaScript(ScriptBase):
    def get_preamble(self):
        return textwrap.dedent(
            """
                ######## snakemake preamble start (automatically inserted, do not edit) ########
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
                ######## snakemake preamble end #########
                """.format(
                JuliaEncoder.encode_namedlist(self.input),
                JuliaEncoder.encode_namedlist(self.output),
                JuliaEncoder.encode_namedlist(self.params),
                JuliaEncoder.encode_namedlist(self.wildcards),
                JuliaEncoder.encode_value(self.threads),
                JuliaEncoder.encode_namedlist(self.log),
                JuliaEncoder.encode_namedlist(
                    {
                        name: value
                        for name, value in self.resources.items()
                        if name != "_cores" and name != "_nodes"
                    }
                ),
                JuliaEncoder.encode_dict(self.config),
                JuliaEncoder.encode_value(self.rulename),
                JuliaEncoder.encode_value(self.bench_iteration),
                JuliaEncoder.encode_value(
                    os.path.dirname(self.path[7:])
                    if self.path.startswith("file://")
                    else os.path.dirname(self.path)
                ),
            ).replace(
                "'", '"'
            )
        )

    def write_script(self, preamble, fd):
        fd.write(preamble.encode())
        fd.write(self.source)

    def execute_script(self, fname, edit=False):
        self._execute_cmd("julia {fname:q}", fname=fname)


def get_source(path, basedir=".", wildcards=None, params=None):
    source = None
    if not path.startswith("http") and not path.startswith("git+file"):
        if path.startswith("file://"):
            path = path[7:]
        elif path.startswith("file:"):
            path = path[5:]
        if not os.path.isabs(path):
            path = smart_join(basedir, path, abspath=True)
        if is_local_file(path):
            path = "file://" + path
    if wildcards is not None and params is not None:
        # Format path if wildcards are given.
        path = format(path, wildcards=wildcards, params=params)
    if path.startswith("file://"):
        sourceurl = "file:" + pathname2url(path[7:])
    elif path.startswith("git+file"):
        source = git_content(path).encode()
        (root_path, file_path, version) = split_git_path(path)
        path = path.rstrip("@" + version)
    else:
        sourceurl = path

    if source is None:
        with urlopen(sourceurl) as source:
            source = source.read()

    language = get_language(path, source)

    return path, source, language


def get_language(path, source):
    import nbformat

    language = None
    if path.endswith(".py"):
        language = "python"
    elif path.endswith(".ipynb"):
        language = "jupyter"
    elif path.endswith(".R"):
        language = "r"
    elif path.endswith(".Rmd"):
        language = "rmarkdown"
    elif path.endswith(".jl"):
        language = "julia"

    # detect kernel language for Jupyter Notebooks
    if language == "jupyter":
        nb = nbformat.reads(source, as_version=nbformat.NO_CONVERT)
        try:
            kernel_language = nb["metadata"]["language_info"]["name"]
        except KeyError as e:
            raise WorkflowError(
                "Notebook metadata is corrupt. Please delete notebook "
                "and recreate it via --edit-notebook."
            )

        language += "_" + kernel_language.lower()

    return language


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
    container_img,
    singularity_args,
    env_modules,
    bench_record,
    jobid,
    bench_iteration,
    cleanup_scripts,
    shadow_dir,
):
    """
    Load a script from the given basedir + path and execute it.
    """
    path, source, language = get_source(path, basedir, wildcards, params)

    exec_class = {
        "python": PythonScript,
        "r": RScript,
        "rmarkdown": RMarkdown,
        "julia": JuliaScript,
    }.get(language, None)
    if exec_class is None:
        raise ValueError(
            "Unsupported script: Expecting either Python (.py), R (.R), RMarkdown (.Rmd) or Julia (.jl) script."
        )

    executor = exec_class(
        path,
        source,
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
        container_img,
        singularity_args,
        env_modules,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
    )
    executor.evaluate()
