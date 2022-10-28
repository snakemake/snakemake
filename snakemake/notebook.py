from abc import abstractmethod
import os, sys
from pathlib import Path
from urllib.error import URLError
import tempfile
import re
import shutil

from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.script import get_source, ScriptBase, PythonScript, RScript
from snakemake.logging import logger
from snakemake.common import is_local_file
from snakemake.common import ON_WINDOWS
from snakemake.sourcecache import SourceCache, infer_source_file
from snakemake.utils import format

KERNEL_STARTED_RE = re.compile(r"Kernel started: (?P<kernel_id>\S+)")
KERNEL_SHUTDOWN_RE = re.compile(r"Kernel shutdown: (?P<kernel_id>\S+)")


class EditMode:
    def __init__(self, server_addr=None, draft_only=False):
        if server_addr is not None:
            self.ip, self.port = server_addr.split(":")
        self.draft_only = draft_only


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
            out.write(nbformat.writes(nb).encode())

    def draft_and_edit(self, listen):
        self.draft()

        self.source = open(self.local_path).read()

        self.evaluate(edit=listen)

    def write_script(self, preamble, fd):
        import nbformat

        nb = nbformat.reads(self.source, as_version=nbformat.NO_CONVERT)

        self.remove_preamble_cell(nb)
        self.insert_preamble_cell(preamble, nb)

        fd.write(nbformat.writes(nb).encode())

    def execute_script(self, fname, edit=None):
        import nbformat

        fname_out = self.log.get("notebook", None)
        if fname_out is None or edit:
            output_parameter = ""
        else:
            fname_out = os.path.join(os.getcwd(), fname_out)
            output_parameter = "--output {fname_out:q}"

        with tempfile.TemporaryDirectory() as tmp:
            if edit is not None:
                assert not edit.draft_only
                logger.info("Opening notebook for editing.")
                cmd = (
                    "jupyter notebook --browser ':' --no-browser --log-level ERROR --ip {edit.ip} --port {edit.port} "
                    "--NotebookApp.quit_button=True {{fname:q}}".format(edit=edit)
                )
            else:
                cmd = (
                    "jupyter-nbconvert --log-level ERROR --execute {output_parameter} "
                    "--to notebook --ExecutePreprocessor.timeout=-1 {{fname:q}}".format(
                        output_parameter=output_parameter,
                    )
                )

            if ON_WINDOWS:
                fname = fname.replace("\\", "/")
                fname_out = fname_out.replace("\\", "/") if fname_out else fname_out

            self._execute_cmd(
                cmd,
                fname_out=fname_out,
                fname=fname,
                additional_envvars={"IPYTHONDIR": tmp},
            )

            if edit:
                logger.info("Saving modified notebook.")
                nb = nbformat.read(fname, as_version=4)

                self.remove_preamble_cell(nb)

                # clean up all outputs
                for cell in nb["cells"]:
                    cell["outputs"] = []

                nbformat.write(nb, self.local_path)

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
    def get_language_name(self):
        ...

    @abstractmethod
    def get_interpreter_exec(self):
        ...


class PythonJupyterNotebook(JupyterNotebook):
    def get_preamble(self):
        preamble_addendum = "import os; os.chdir(r'{cwd}');".format(cwd=os.getcwd())

        return PythonScript.generate_preamble(
            self.path,
            self.cache_path,
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

    def get_language_name(self):
        return "python"

    def get_interpreter_exec(self):
        return "python"


class RJupyterNotebook(JupyterNotebook):
    def get_preamble(self):
        preamble_addendum = "setwd('{cwd}');".format(cwd=os.getcwd())

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
            preamble_addendum=preamble_addendum,
        )

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
    edit,
    runtime_sourcecache_path,
):
    """
    Load a script from the given basedir + path and execute it.
    """
    draft = False
    path = format(path, wildcards=wildcards, params=params)
    if edit is not None:
        if is_local_file(path):
            if not os.path.isabs(path):
                local_path = os.path.join(basedir, path)
            else:
                local_path = path
            if not os.path.exists(local_path):
                # draft the notebook, it does not exist yet
                language = None
                draft = True
                path = "file://{}".format(os.path.abspath(local_path))
                if path.endswith(".py.ipynb"):
                    language = "jupyter_python"
                elif path.endswith(".r.ipynb"):
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
            path, SourceCache(runtime_sourcecache_path), basedir, wildcards, params
        )
    else:
        source = None
        cache_path = None
        is_local = True
        path = infer_source_file(path)

    exec_class = get_exec_class(language)

    executor = exec_class(
        path,
        cache_path,
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

    if edit is None:
        executor.evaluate(edit=edit)
    elif edit.draft_only:
        executor.draft()
        msg = "Generated skeleton notebook:\n{} ".format(path)
        if conda_env and not container_img:
            msg += (
                "\n\nEditing with VSCode:\nOpen notebook, run command 'Select notebook kernel' (Ctrl+Shift+P or Cmd+Shift+P), and choose:"
                "\n{}\n".format(
                    str(Path(conda_env) / "bin" / executor.get_interpreter_exec())
                )
            )
            msg += (
                "\nEditing with Jupyter CLI:"
                "\nconda activate {}\njupyter notebook {}\n".format(conda_env, path)
            )
        logger.info(msg)
    elif draft:
        executor.draft_and_edit(listen=edit)
    else:
        executor.evaluate(edit=edit)
