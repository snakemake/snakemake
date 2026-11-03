from typing import Optional
from snakemake.settings.types import NotebookEditMode
from snakemake.executors.local import RunArgs
from typing import Dict
from abc import abstractmethod
import os
from pathlib import Path
import shutil
import tempfile
import re
from textwrap import dedent

from snakemake.exceptions import WorkflowError
from snakemake.script import get_source, ScriptBase, PythonScript, RScript
from snakemake.logging import logger
from snakemake.common.misc import is_local_file
from snakemake.common.constants import ON_WINDOWS
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
                has_papermill = (
                    self.run_args.software_env
                    and self.run_args.software_env.contains_executable("papermill")
                ) or shutil.which("papermill") is not None
                if has_papermill:
                    if fname_out is None:
                        output_parameter = fname
                    else:
                        output_parameter = "{fname_out}"
                    cmd = "papermill --log-level ERROR {fname:q} " f"{output_parameter}"
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


class MarimoNotebook(PythonScript):
    editable = True

    def get_interpreter_exec(self):
        return "marimo"

    def insert_preamble_cell(self, preamble, notebook):
        fixed_preamble = re.sub(
            "__real_file__ = __file__", "__real_file__ = __name__", preamble
        ).replace("\\", r"\\")

        preamble_cell = "\n".join(
            [
                "@app.cell(hide_code=True)",
                "def snakemake_preamble():",
                *[f"    {line}" for line in fixed_preamble.splitlines()],
                "    return\n\n",
            ]
        )

        notebook_with_preamble = re.sub(
            '(if __name__ == "__main__":)', rf"{preamble_cell}\g<1>", notebook
        )

        if not re.search("def snakemake_preamble", notebook_with_preamble):
            raise ValueError(dedent("""\
                    Could not inject Snakemake preamble into marimo notebook.
                    Are you sure it is a valid marimo notebook?
                    Hint: Use `marimo check` to check the format of the file.
                    """))

        return notebook_with_preamble

    def remove_preamble_cell(self, notebook):
        return re.sub(
            r"""
            ^@app\.cell.*\n
            def\ snakemake_preamble.*\n
            (?:(?!\s{4}return).*\n)+
            \s{4}return.*$
            \s+
            """,
            "",
            notebook,
            flags=re.M | re.X,
        )

    def draft(self):
        minimal_notebook = dedent("""\
            import marimo

            app = marimo.App()

            @app.cell
            def _():
                # start coding here
                return

            if __name__ == "__main__":
                app.run()
        """)
        notebook_with_preamble = self.insert_preamble_cell(
            self.get_preamble(), minimal_notebook
        )

        os.makedirs(os.path.dirname(self.local_path), exist_ok=True)

        with open(self.local_path, "wb") as out:
            out.write(notebook_with_preamble.encode())

    def draft_and_edit(self, listen):
        self.draft()
        self.source = open(self.local_path).read()
        self.evaluate(edit=listen)

    def write_script(self, preamble, fd):
        source_without_preamble = self.remove_preamble_cell(self.source)
        source_with_preamble = self.insert_preamble_cell(
            preamble, source_without_preamble
        )
        fd.write(source_with_preamble.encode())

    def execute_script(self, fname, edit=None):
        fname_out = self.log.get("notebook", None)
        if fname_out is not None:
            if not fname_out.endswith(".html"):
                raise ValueError(
                    "Only `html` is supported for saving processed marimo notebooks."
                )
            fname_out = os.path.abspath(fname_out)
        export_cmd = "marimo export html {fname:q} -o {fname_out:q}"

        if edit is not None:
            assert not edit.draft_only
            cmd = "marimo edit {fname:q} --port {edit.port} --host {edit.ip}"
            if fname_out is not None:
                cmd += f" && {export_cmd}"
        elif fname_out is not None:
            cmd = export_cmd
        else:
            cmd = "python {fname:q}"

        self._execute_cmd(cmd, fname=fname, fname_out=fname_out, edit=edit)

        if edit:
            logger.info("Saving modified notebook.")
            with open(fname, "r") as file:
                notebook = self.remove_preamble_cell(file.read())
            with open(self.local_path, "w") as out_file:
                out_file.write(notebook)


def get_exec_class(language):
    exec_class = {
        "jupyter_python": PythonJupyterNotebook,
        "jupyter_r": RJupyterNotebook,
        "marimo_python": MarimoNotebook,
    }.get(language, None)
    if exec_class is None:
        raise ValueError(
            "Unsupported notebook: Expecting Jupyter (.ipynb) or marimo (.py) notebook."
        )
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
                local_path = Path(
                    run_args.basedir.join(path).get_path_or_uri(secret_free=True)
                )
            else:
                local_path = Path(path)
            if not local_path.exists():
                # draft the notebook, it does not exist yet
                language = None
                draft = True
                path = f"file://{local_path.absolute()}"
                suffixes = local_path.suffixes
                if suffixes[-2:] == [".py", ".ipynb"]:
                    language = "jupyter_python"
                elif suffixes[-2:] == [".r", ".ipynb"]:
                    language = "jupyter_r"
                elif path.endswith(".marimo.py"):
                    language = "marimo_python"
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
