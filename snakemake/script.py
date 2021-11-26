__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import inspect
import itertools
import os
from snakemake import sourcecache
from snakemake.sourcecache import (
    LocalSourceFile,
    SourceCache,
    SourceFile,
    infer_source_file,
)
import tempfile
import textwrap
import sys
import pickle
import subprocess
import collections
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Tuple, Pattern, Union, Optional
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
PathLike = Union[str, Path, os.PathLike]


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
        return _log_shell_redirect(self.log, stdout, stderr, append)


def _log_shell_redirect(
    log: Optional[PathLike],
    stdout: bool = True,
    stderr: bool = True,
    append: bool = False,
) -> str:
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
    if not log:
        return ""
    lookup = {
        (True, True, True): " >> {0} 2>&1",
        (True, False, True): " >> {0}",
        (False, True, True): " 2>> {0}",
        (True, True, False): " > {0} 2>&1",
        (True, False, False): " > {0}",
        (False, True, False): " 2> {0}",
    }
    return lookup[(stdout, stderr, append)].format(str(log))


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
        elif isinstance(value, Path):
            return repr(str(value))
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
        elif isinstance(value, Path):
            return repr(str(value))
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
        conda_base_path,
        container_img,
        singularity_args,
        env_modules,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
        is_local,
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
        self.conda_base_path = conda_base_path
        self.container_img = container_img
        self.singularity_args = singularity_args
        self.env_modules = env_modules
        self.bench_record = bench_record
        self.jobid = jobid
        self.bench_iteration = bench_iteration
        self.cleanup_scripts = cleanup_scripts
        self.shadow_dir = shadow_dir
        self.is_local = is_local

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
                suffix="." + self.path.get_filename(), dir=dir_, delete=False
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
        assert self.is_local
        return self.path.get_path_or_uri()

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
            conda_base_path=self.conda_base_path,
            container_img=self.container_img,
            shadow_dir=self.shadow_dir,
            env_modules=self.env_modules,
            singularity_args=self.singularity_args,
            resources=self.resources,
            threads=self.threads,
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
        is_local,
        preamble_addendum="",
    ):
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
            path.get_basedir().get_path_or_uri(),
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
        if is_local:
            searchpath += ", " + repr(path.get_basedir().get_path_or_uri())

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

        if isinstance(self.path, LocalSourceFile):
            file_override = os.path.realpath(self.path.get_path_or_uri())
        else:
            file_override = self.path.get_path_or_uri()
        preamble_addendum = (
            "__real_file__ = __file__; __file__ = {file_override};".format(
                file_override=repr(file_override)
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
            self.is_local,
            preamble_addendum=preamble_addendum,
        )

    def write_script(self, preamble, fd):
        fd.write(preamble.encode())
        fd.write(self.source.encode())

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
        self._execute_cmd(
            "{py_exec} {fname:q}", py_exec=py_exec, fname=fname, is_python_script=True
        )


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
            REncoder.encode_value(path.get_basedir().get_path_or_uri()),
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
        fd.write(self.source.encode())

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
            REncoder.encode_value(self.path.get_basedir().get_path_or_uri()),
        )

    def write_script(self, preamble, fd):
        # Insert Snakemake object after the RMarkdown header
        code = self.source
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
        fd.write(code[pos:].encode())

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
                JuliaEncoder.encode_value(self.path.get_basedir().get_path_or_uri()),
            ).replace(
                "'", '"'
            )
        )

    def write_script(self, preamble, fd):
        fd.write(preamble.encode())
        fd.write(self.source.encode())

    def execute_script(self, fname, edit=False):
        self._execute_cmd("julia {fname:q}", fname=fname)


class RustScript(ScriptBase):
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
        is_local,
        preamble_addendum="",
    ):
        # snakemake's namedlists will be encoded as a dict
        # which stores the not-named items at the key "positional"
        # and unpacks named items into the dict
        def encode_namedlist(values):
            values = list(values)
            if len(values) == 0:
                return dict(positional=[])
            positional = [val for key, val in values if not key]
            return dict(
                positional=positional, **{key: val for key, val in values if key}
            )

        snakemake = dict(
            input=encode_namedlist(input_._plainstrings()._allitems()),
            output=encode_namedlist(output._plainstrings()._allitems()),
            params=encode_namedlist(params.items()),
            wildcards=encode_namedlist(wildcards.items()),
            threads=threads,
            resources=encode_namedlist(
                {
                    name: value
                    for (name, value) in resources.items()
                    if name != "_cores" and name != "_nodes"
                }.items()
            ),
            log=encode_namedlist(log._plainstrings()._allitems()),
            config=encode_namedlist(config.items()),
            rulename=rulename,
            bench_iteration=bench_iteration,
            scriptdir=path.get_basedir().get_path_or_uri(),
        )

        import json

        json_string = json.dumps(dict(snakemake))

        # Obtain search path for current snakemake module.
        # We append it at the end (as a fallback).
        searchpath = SNAKEMAKE_SEARCHPATH
        if container_img is not None:
            searchpath = singularity.SNAKEMAKE_MOUNTPOINT
        searchpath = repr(searchpath)
        # For local scripts, add their location to the path in case they use path-based imports
        if is_local:
            searchpath += ", " + repr(path.get_basedir().get_path_or_uri())

        return textwrap.dedent(
            """
            json_typegen::json_typegen!("Snakemake", r###"{json_string}"###, {{
                "/bench_iteration": {{
                   "use_type": "Option<usize>"
                }},
                "/input/positional": {{
                    "use_type": "Vec<String>"
                }},
                "/output/positional": {{
                    "use_type": "Vec<String>"
                }},
                "/log/positional": {{
                    "use_type": "Vec<String>"
                }},
                "/wildcards/positional": {{
                    "use_type": "Vec<String>"
                }},
            }});
            
            pub struct Iter<'a, T>(std::slice::Iter<'a, T>);
            impl<'a, T> Iterator for Iter<'a, T> {{
                type Item = &'a T;
                
                fn next(&mut self) -> Option<Self::Item> {{
                    self.0.next()
                }}
            }}
            macro_rules! impl_iter {{
                ($($s:ty),+) => {{
                    $(
                        impl IntoIterator for $s {{
                            type Item = String;
                            type IntoIter = std::vec::IntoIter<Self::Item>;
            
                            fn into_iter(self) -> Self::IntoIter {{
                                self.positional.into_iter()
                            }}
                        }}
            
                        impl<'a> IntoIterator for &'a $s {{
                            type Item = &'a String;
                            type IntoIter = Iter<'a, String>;
            
                            fn into_iter(self) -> Self::IntoIter {{
                                Iter(self.positional.as_slice().into_iter())
                            }}
                        }}
                    )+
                }};
            }}
            
            macro_rules! impl_index {{
                ($($s:ty),+) => {{
                    $(
                    impl std::ops::Index<usize> for $s {{
                        type Output = String;
            
                        fn index(&self, index: usize) -> &Self::Output {{
                            &self.positional[index]
                        }}
                    }}
                    )+
                }}
            }}
            
            
            impl_iter!(Input, Output, Wildcards, Log);
            impl_index!(Input, Output, Wildcards, Log);
            
            impl Snakemake {{
                #[allow(dead_code)]
                fn redirect_stderr<P: AsRef<std::path::Path>>(
                    &self,
                    path: P,
                ) -> anyhow::Result<gag::Redirect<std::fs::File>> {{
                    let log = std::fs::OpenOptions::new()
                        .truncate(true)
                        .read(true)
                        .create(true)
                        .write(true)
                        .open(path)?;
                    Ok(gag::Redirect::stderr(log)?)
                }}
                
                #[allow(dead_code)]
                fn redirect_stdout<P: AsRef<std::path::Path>>(
                    &self,
                    path: P,
                ) -> anyhow::Result<gag::Redirect<std::fs::File>> {{
                    let log = std::fs::OpenOptions::new()
                        .truncate(true)
                        .read(true)
                        .create(true)
                        .write(true)
                        .open(path)?;
                    Ok(gag::Redirect::stdout(log)?)
                }}
                
                fn setup_path(&self) -> anyhow::Result<()> {{
                    use std::env;
                    if let Some(path) = env::var_os("PATH") {{
                        let mut paths = env::split_paths(&path).collect::<Vec<_>>();
                        paths.push(std::path::PathBuf::from("{searchpath}"));
                        let new_path = env::join_paths(paths)?;
                        env::set_var("PATH", &new_path);
                    }}
                    Ok(())
                }}
            }}
            
            lazy_static::lazy_static! {{
                // https://github.com/rust-lang-nursery/lazy-static.rs/issues/153
                #[allow(non_upper_case_globals)]
                static ref snakemake: Snakemake = {{
                    let s: Snakemake = serde_json::from_str(r###"{json_string}"###).expect("Failed parsing snakemake JSON");
                    s.setup_path().expect("Failed setting PATH");
                    s
                }};
            }}
            // TODO include addendum, if any {{preamble_addendum}}
            """
        ).format(
            searchpath=searchpath,
            json_string=json_string,
            preamble_addendum=preamble_addendum,
        )

    def get_preamble(self):
        preamble_addendum = ""

        preamble = RustScript.generate_preamble(
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
            self.is_local,
            preamble_addendum=preamble_addendum,
        )
        return preamble

    def write_script(self, preamble, fd):
        content = self.combine_preamble_and_source(preamble)
        fd.write(content.encode())

    def execute_script(self, fname, edit=False):
        deps = self.default_dependencies()
        ftrs = self.default_features()
        self._execute_cmd(
            "rust-script -d {deps} --features {ftrs} {fname:q} ",
            fname=fname,
            deps=deps,
            ftrs=ftrs,
        )

    def combine_preamble_and_source(self, preamble: str) -> str:
        """The manifest info needs to be moved to before the preamble.
        Also, because rust-scipt relies on inner docs, there can't be an empty line
        between the manifest and preamble.
        """
        manifest, src = RustScript.extract_manifest(self.source)
        return manifest + preamble.lstrip("\r\n") + src

    @staticmethod
    def default_dependencies() -> str:
        return " -d ".join(
            [
                "anyhow=1",
                "serde_json=1",
                "serde=1",
                "serde_derive=1",
                "lazy_static=1.4",
                "json_typegen=0.6",
                "gag=1",
            ]
        )

    @staticmethod
    def default_features() -> str:
        return ",".join(["serde/derive"])

    @staticmethod
    def extract_manifest(source: str) -> Tuple[str, str]:
        # we have no need for the shebang for now given the way we run the script
        _, src = RustScript._strip_shebang(source)
        manifest, src = RustScript._strip_manifest(src)

        return manifest, src

    @staticmethod
    def _strip_shebang(src: str) -> Tuple[str, str]:
        """From https://github.com/fornwall/rust-script/blob/ce508bad02a11d574657d2f1debf7e73fca2bf6e/src/manifest.rs#L312-L320"""
        rgx = re.compile(r"^#![^\[].*?(\r\n|\n)")
        return strip_re(rgx, src)

    @staticmethod
    def _strip_manifest(src: str) -> Tuple[str, str]:
        """From https://github.com/fornwall/rust-script/blob/ce508bad02a11d574657d2f1debf7e73fca2bf6e/src/manifest.rs#L405-L411"""
        manifest, remainder = RustScript._strip_single_line_manifest(src)
        if not manifest:
            manifest, remainder = RustScript._strip_code_block_manifest(src)
        return manifest, remainder

    @staticmethod
    def _strip_single_line_manifest(src: str) -> Tuple[str, str]:
        """From https://github.com/fornwall/rust-script/blob/ce508bad02a11d574657d2f1debf7e73fca2bf6e/src/manifest.rs#L618-L632"""
        rgx = re.compile(r"^\s*//\s*cargo-deps\s*:(.*?)(\r\n|\n)", flags=re.IGNORECASE)
        return strip_re(rgx, src)

    @staticmethod
    def _strip_code_block_manifest(src: str) -> Tuple[str, str]:
        """From https://github.com/fornwall/rust-script/blob/ce508bad02a11d574657d2f1debf7e73fca2bf6e/src/manifest.rs#L634-L664
        We need to find the first `/*!` or `//!` that *isn't* preceeded by something
        that would make it apply to anything other than the crate itself. Because we
        can't do this accurately, we'll just require that the doc comment is the
        *first* thing in the file (after the optional shebang, which should already
        have been stripped).
        """
        crate_comment_re = re.compile(
            r"^\s*(/\*!|//([!/]))(.*?)(\r\n|\n)", flags=re.MULTILINE
        )
        # does src start with a crate comment?
        match = crate_comment_re.match(src)
        if not match:
            return "", src
        end_of_comment = match.end()
        # find end of crate comment
        while match is not None:
            end_of_comment = match.end()
            match = crate_comment_re.match(src, pos=end_of_comment)

        crate_comment = src[:end_of_comment]
        found_code_block_open = False
        code_block_open_re = re.compile(r"```\s*cargo")
        found_code_block_close = False
        code_block_close_re = re.compile(r"```")
        for line in crate_comment.splitlines():
            if not found_code_block_open:
                m = code_block_open_re.search(line)
                if m:
                    found_code_block_open = True
            else:
                m = code_block_close_re.search(line)
                if m:
                    found_code_block_close = True
                    break

        crate_comment_has_manifest = found_code_block_open and found_code_block_close
        if crate_comment_has_manifest:
            return crate_comment, src[end_of_comment:]
        else:
            return "", src


def strip_re(regex: Pattern, s: str) -> Tuple[str, str]:
    """Strip a substring matching a regex from a string and return the stripped part
    and the remainder of the original string.
    Returns an empty string and the original string if the regex is not found
    """
    rgx = re.compile(regex)
    match = rgx.search(s)
    if match:
        head, tail = s[: match.end()], s[match.end() :]
    else:
        head, tail = "", s

    return head, tail


def get_source(
    path,
    sourcecache: sourcecache.SourceCache,
    basedir=None,
    wildcards=None,
    params=None,
):
    if wildcards is not None and params is not None:
        if isinstance(path, SourceFile):
            path = path.get_path_or_uri()
        # Format path if wildcards are given.
        path = infer_source_file(format(path, wildcards=wildcards, params=params))

    if basedir is not None:
        basedir = infer_source_file(basedir)

    source_file = infer_source_file(path, basedir)
    with sourcecache.open(source_file) as f:
        source = f.read()

    language = get_language(source_file, source)

    is_local = isinstance(source_file, LocalSourceFile)

    return source_file, source, language, is_local


def get_language(source_file, source):
    import nbformat

    filename = source_file.get_filename()

    language = None
    if filename.endswith(".py"):
        language = "python"
    elif filename.endswith(".ipynb"):
        language = "jupyter"
    elif filename.endswith(".R"):
        language = "r"
    elif filename.endswith(".Rmd"):
        language = "rmarkdown"
    elif filename.endswith(".jl"):
        language = "julia"
    elif filename.endswith(".rs"):
        language = "rust"

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
    conda_base_path,
    container_img,
    singularity_args,
    env_modules,
    bench_record,
    jobid,
    bench_iteration,
    cleanup_scripts,
    shadow_dir,
    runtime_sourcecache_path,
):
    """
    Load a script from the given basedir + path and execute it.
    """

    path, source, language, is_local = get_source(
        path, SourceCache(runtime_sourcecache_path), basedir, wildcards, params
    )

    exec_class = {
        "python": PythonScript,
        "r": RScript,
        "rmarkdown": RMarkdown,
        "julia": JuliaScript,
        "rust": RustScript,
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
        conda_base_path,
        container_img,
        singularity_args,
        env_modules,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
        is_local,
    )
    executor.evaluate()
