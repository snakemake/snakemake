import os
import pickle
import tempfile
import textwrap
from urllib.error import URLError
from abc import ABC, abstractmethod

import nbformat

from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.common import escape_backslash, SNAKEMAKE_SEARCHPATH
from snakemake import singularity


class ScriptBase(ABC):
    def __init__(
        self,
        path, source,
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
        singularity_img,
        singularity_args,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir
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
        self.singularity_img = singularity_img
        self.singularity_args = singularity_args
        self.bench_record = bench_record
        self.jobid = jobid
        self.bench_iteration = bench_iteration
        self.cleanup_scripts = cleanup_scripts
        self.shadow_dir = shadow_dir

    def evaluate(self):
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
            self.execute_script(fd.name)
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

    @abstractmethod
    def get_preamble(self):
        ...

    @abstractmethod
    def write_script(self, preamble, fd):
        ...

    @abstractmethod
    def execute_script(self, fname):
        ...


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


class JupyterNotebook(ScriptBase):
    def _get_kernel_language(self):
        nb = nbformat.reads(self.source, as_version=nbformat.NO_CONVERT)
        return nb['metadata']['language_info']['name']

    def _get_python_preamble(self):
        wrapper_path = self.path[7:] if self.path.startswith("file://") else self.path
        snakemake = Snakemake(
            self.input,
            self.output,
            self.params,
            self.wildcards,
            self.threads,
            self.resources,
            self.log,
            self.config,
            self.rulename,
            self.bench_iteration,
            os.path.dirname(wrapper_path),
        )
        snakemake = pickle.dumps(snakemake)
        # Obtain search path for current snakemake module.
        # The module is needed for unpickling in the script.
        # We append it at the end (as a fallback).
        searchpath = SNAKEMAKE_SEARCHPATH
        if self.singularity_img is not None:
            searchpath = singularity.SNAKEMAKE_MOUNTPOINT
        searchpath = '"{}"'.format(searchpath)
        # For local scripts, add their location to the path in case they use path-based imports
        if self.path.startswith("file://"):
            searchpath += ', "{}"'.format(os.path.dirname(self.path[7:]))

        # nbconvert sets cwd to notebook directory.
        # This is problematic because we create a temporary file.
        preamble_addendum = "import os; os.chdir('{cwd}');".format(
            cwd=os.getcwd()
        )

        return textwrap.dedent(
            """
        ######## Snakemake header ########
        import sys; sys.path.extend([{searchpath}]); import pickle; snakemake = pickle.loads({snakemake}); from snakemake.logging import logger; logger.printshellcmds = {printshellcmds}; {preamble_addendum}
        ######## Original script #########
        """
        ).format(
            searchpath=escape_backslash(searchpath),
            snakemake=snakemake,
            printshellcmds=logger.printshellcmds,
            preamble_addendum=preamble_addendum,
        )

    def _get_r_preamble(self):
        # nbconvert sets cwd to notebook directory.
        # This is problematic because we create a temporary file.
        preamble_addendum = "setwd('{cwd}');".format(
            cwd=os.getcwd()
        )

        return textwrap.dedent(
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
        {preamble_addendum}

        ######## Original script #########
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
            preamble_addendum=preamble_addendum
        )

    def get_preamble(self):
        kernel_language = self._get_kernel_language()
        if kernel_language == 'python':
            return self._get_python_preamble()
        elif kernel_language == 'R':
            return self._get_r_preamble()
        else:
            raise ValueError(
                "Unsupported Jupyter kernel: Expecting one of the following: Python"
            )

    def write_script(self, preamble, fd):
        nb = nbformat.reads(self.source, as_version=4)  # nbformat.NO_CONVERT

        preamble_cell = nbformat.v4.new_code_cell(preamble)
        nb["cells"].insert(0, preamble_cell)

        fd.write(nbformat.writes(nb).encode())

    def execute_script(self, fname):
        # execute notebook
        tmp_output = "{fname}.processed.ipynb".format(fname=fname)
        shell(
            "jupyter nbconvert --execute --output {tmp_output:q} --to notebook --ExecutePreprocessor.timeout=-1 {fname:q}",
            bench_record=self.bench_record,
            tmp_output=tmp_output,
        )

        # determine whether to save output
        notebook_relpath = self.output.get("notebook_output", None)

        if notebook_relpath is not None:
            # determine output format
            _, ext = os.path.splitext(notebook_relpath)
            output_format = {
                ".ipynb": "notebook",
                ".html": "html",
                ".tex": "latex",
                ".pdf": "pdf",
                # ".slides": "slides",
                ".md": "markdown",
                ".txt": "asciidoc",
                ".rst": "rst",
                # ".py": "script",
            }.get(ext, None)

            if output_format is None:
                raise WorkflowError(
                    "Invalid Jupyter Notebook output format: '{ext}'".format(
                        ext=ext
                    )
                )

            # remove preamble from output
            nb = nbformat.read(tmp_output, as_version=nbformat.NO_CONVERT)

            nb["cells"].pop(0)

            with open(tmp_output, "w") as fd:
                nbformat.write(nb, fd)

            # save to destination
            notebook_output = os.path.join(os.getcwd(), notebook_relpath)

            if output_format == "notebook":
                os.rename(tmp_output, notebook_output)
            else:
                # convert to required format
                shell(
                    "jupyter nbconvert --output {notebook_output:q} --to {output_format:q} --ExecutePreprocessor.timeout=-1 {tmp_output:q}",
                    notebook_output=notebook_output,
                    output_format=output_format,
                    tmp_output=tmp_output,
                )


def get_executor_class(language):
    ExecClass = {
        'jupyter': JupyterNotebook
    }.get(language, None)

    if ExecClass is None:
        raise ValueError(
            "Unsupported script: Expecting either Python (.py), Jupyter Notebook (.ipynb), R (.R), RMarkdown (.Rmd) or Julia (.jl) script."
        )
    return ExecClass
