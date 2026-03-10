from typing import Optional
from snakemake.settings.types import NotebookEditMode
from snakemake.executors.local import RunArgs
from typing import Dict
from abc import abstractmethod
import os
from pathlib import Path
import subprocess as sp
import shutil
import tempfile
import re

from snakemake.exceptions import WorkflowError
from snakemake.script import get_source, ScriptBase, PythonScript, RScript
from snakemake.logging import logger
from snakemake.common import is_local_file
from snakemake.common import ON_WINDOWS
from snakemake.sourcecache import SourceCache, infer_source_file
from snakemake.utils import format

KERNEL_STARTED_RE = re.compile(r"Kernel started: (?P<kernel_id>\S+)")
KERNEL_SHUTDOWN_RE = re.compile(r"Kernel shutdown: (?P<kernel_id>\S+)")

NBFORMAT_VERSION = 4


def get_cell_sources(source):
    import nbformat

    nb = nbformat.reads(source, as_version=nbformat.NO_CONVERT)

    return [cell["source"] for cell in nb["cells"]]


class JupyterNotebook(ScriptBase):
    editable = True

    def draft(self):
        import nbformat

        preamble = self.get_preamble()
        nb = nbformat.v4.new_notebook()
        self.insert_preamble_cell(preamble, nb)

        nb["cells"].append(nbformat.v4.new_code_cell("# start coding here"))
        nb["metadata"] = {"language_info": {"name": self.get_language_name()}}

        os.makedirs(os.path.dirname(self.local_path), exist_ok=True)

        with open(self.local_path, "wb") as out:
            out.write(nbformat.writes(nb, version=NBFORMAT_VERSION).encode())

    def draft_and_edit(self, listen):
        self.draft()

        self.source = open(self.local_path).read()

        self.evaluate(edit=listen)

    def write_script(self, preamble, fd):
        import nbformat

        nb = nbformat.reads(self.source, as_version=NBFORMAT_VERSION)

        self.remove_preamble_cell(nb)
        self.insert_preamble_cell(preamble, nb)

        fd.write(nbformat.writes(nb, version=NBFORMAT_VERSION).encode())

    def execute_script(self, fname, edit: Optional[NotebookEditMode] = None):
        import nbformat

        fname_out = self.run_args.log.get("notebook", None)

        with tempfile.TemporaryDirectory() as tmp:
            if edit is not None:
                assert not edit.draft_only
                logger.info(f"Opening notebook for editing at {edit.ip}:{edit.port}")
                cmd = (
                    f"jupyter notebook --browser ':' --no-browser --log-level ERROR --ip {edit.ip} --port {edit.port} "
                    "--ServerApp.quit_button=True {fname:q}"
                )
            else:
                has_papermill = (self.run_args.software_env and self.run_args.software_env.contains_executable("papermill")) or shutil.which("papermill") is not None
                if has_papermill:
                    if fname_out is None:
                        output_parameter = fname
                    else:
                        output_parameter = "{fname_out}"
                    cmd = (
                        "papermill --log-level ERROR {fname:q} "
                        f"{output_parameter}"
                    )
                else:
                    if fname_out is None:
                        output_parameter = f"--output '{tmp}/notebook.ipynb'"
                    else:
                        fname_out = os.path.abspath(fname_out)
                        output_parameter = "--output {fname_out:q}"

                    cmd = (
                        f"jupyter-nbconvert --log-level ERROR --execute {output_parameter} "
                        "--to notebook --ExecutePreprocessor.timeout=-1 {fname:q}"
                    )

            if ON_WINDOWS:
                fname = fname.replace("\\", "/")
                fname_out = fname_out.replace("\\", "/") if fname_out else fname_out

            self._execute_cmd(
                cmd,
                fname_out=fname_out,
                fname=fname,
                additional_envvars={"IPYTHONDIR": tmp},
                is_python_script=True,
            )

            if edit:
                if fname_out is not None:
                    # store log file (executed notebook) in requested path
                    shutil.copyfile(fname, fname_out)

                logger.info("Saving modified notebook.")
                nb = nbformat.read(fname, as_version=NBFORMAT_VERSION)

                self.remove_preamble_cell(nb)

                # clean up all outputs
                for cell in nb["cells"]:
                    if "outputs" in cell:
                        cell["outputs"] = []
                    if "execution_count" in cell:
                        cell["execution_count"] = None

                nbformat.write(nb, self.local_path, version=NBFORMAT_VERSION)

    def insert_preamble_cell(self, preamble, notebook):
        import nbformat

        preamble_cell = nbformat.v4.new_code_cell(preamble)
        preamble_cell["metadata"]["tags"] = ["snakemake-job-properties"]
        notebook["cells"].insert(0, preamble_cell)

    def remove_preamble_cell(self, notebook):
        preambles = [
            i
            for i, cell in enumerate(notebook["cells"])
            if "snakemake-job-properties" in cell["metadata"].get("tags", [])
        ]
        if len(preambles) > 1:
            raise WorkflowError(
                "More than one snakemake preamble cell found in notebook. "
                "Please clean up the notebook first, by removing all or all but one of them."
            )
        elif len(preambles) == 1:
            preamble = preambles[0]
            # remove old preamble
            del notebook["cells"][preamble]

    @abstractmethod
    def get_language_name(self): ...

    @abstractmethod
    def get_interpreter_exec(self): ...


class PythonJupyterNotebook(JupyterNotebook, PythonScript):

    def preamble_addendum(self):
        return f"import os; os.chdir(r'{os.getcwd()}');"

    def get_language_name(self):
        return "python"

    def get_interpreter_exec(self):
        return "python"


class RJupyterNotebook(JupyterNotebook, RScript):
    def preamble_addendum(self) -> str:
        return f"setwd('{os.getcwd()}');"

    def get_language_name(self):
        return "r"

    def get_interpreter_exec(self):
        return "RScript"


def get_exec_class(language):
    exec_class = {
        "jupyter_python": PythonJupyterNotebook,
        "jupyter_r": RJupyterNotebook,
    }.get(language, None)
    if exec_class is None:
        raise ValueError("Unsupported notebook: Expecting Jupyter Notebook (.ipynb).")
    return exec_class


def notebook(
    path,
    run_args: RunArgs,
    config: Dict,
) -> None:
    """
    Load a notebook from the given basedir + path and execute it.
    """
    draft = False
    if isinstance(path, Path):
        path = str(path)
    path = format(path, wildcards=run_args.wildcards, params=run_args.params)
    if run_args.edit_notebook is not None:
        if is_local_file(path):
            if not os.path.isabs(path):
                local_path = Path(run_args.basedir.join(path).get_path_or_uri(secret_free=True))
            else:
                local_path = Path(path)
            if not local_path.exists():
                # draft the notebook, it does not exist yet
                language = None
                draft = True
                path = f"file://{local_path.absolute()}"
                suffixes = path.suffixes
                if suffixes[-2:] == [".py", ".ipynb"]:
                    language = "jupyter_python"
                elif suffixes[-2:] == [".r", ".ipynb"]:
                    language = "jupyter_r"
                else:
                    raise WorkflowError(
                        "Notebook to edit has to end on .py.ipynb or .r.ipynb in order "
                        "to decide which programming language shall be used."
                    )
        else:
            raise WorkflowError(
                "Notebook {} is not local, but edit mode is only allowed for "
                "local notebooks.".format(path)
            )

    if not draft:
        path, source, language, is_local, cache_path = get_source(
            path,
            SourceCache(run_args.cache_path, run_args.runtime_cache_path),
            run_args.basedir,
            run_args.wildcards,
            run_args.params,
        )
    else:
        source = None
        cache_path = None
        is_local = True
        path = infer_source_file(path)

    exec_class = get_exec_class(language)
    if exec_class is None:
        raise ValueError(
            "Unsupported notebook: Has to be either an R (.r.ipynb) or Python (.py.ipynb) notebook."
        )

    executor = exec_class(
        path,
        cache_path=cache_path,
        source=source,
        is_local=is_local,
        run_args=run_args,
        config=config,
    )

    if run_args.edit_notebook is None:
        executor.evaluate(edit=run_args.edit_notebook)
    elif run_args.edit_notebook.draft_only:
        executor.draft()
        msg = f"Generated skeleton notebook:\n{path} "
        # TODO provide hints how to start notebook, also within the annotated
        # software envs.
        logger.info(msg)
    elif draft:
        executor.draft_and_edit(listen=run_args.edit_notebook)
    else:
        executor.evaluate(edit=run_args.edit_notebook)
