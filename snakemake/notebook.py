import os, sys
from urllib.error import URLError
from urllib.parse import urlparse
import tempfile
import re
import shutil

from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.script import get_source, ScriptBase, PythonScript, RScript
from snakemake.logging import logger

KERNEL_STARTED_RE = re.compile(r"Kernel started: (?P<kernel_id>\S+)")
KERNEL_SHUTDOWN_RE = re.compile(r"Kernel shutdown: (?P<kernel_id>\S+)")


class Listen:
    def __init__(self, arg):
        self.ip, self.port = arg.split(":")


class JupyterNotebook(ScriptBase):

    editable = True

    def draft(self, listen):
        import nbformat

        preamble = self.get_preamble()
        nb = nbformat.v4.new_notebook()
        self.insert_preamble_cell(preamble, nb)

        nb["cells"].append(nbformat.v4.new_code_cell("# start coding here"))

        os.makedirs(os.path.dirname(self.local_path), exist_ok=True)

        with open(self.local_path, "wb") as out:
            out.write(nbformat.writes(nb).encode())

        self.source = open(self.local_path).read()

        self.evaluate(edit=listen)

    def write_script(self, preamble, fd):
        import nbformat

        nb = nbformat.reads(self.source, as_version=4)  # nbformat.NO_CONVERT

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

        if edit is not None:
            logger.info("Opening notebook for editing.")
            cmd = (
                "jupyter notebook --browser ':' --no-browser --log-level ERROR --ip {edit.ip} --port {edit.port} "
                "--NotebookApp.quit_button=True {{fname:q}}".format(edit=edit)
            )
        else:
            cmd = (
                "jupyter-nbconvert --log-level ERROR --execute {output_parameter} "
                "--to notebook --ExecutePreprocessor.timeout=-1 {{fname:q}}".format(
                    output_parameter=output_parameter
                )
            )

        self._execute_cmd(cmd, fname_out=fname_out, fname=fname)

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


class PythonJupyterNotebook(JupyterNotebook):
    def get_preamble(self):
        preamble_addendum = "import os; os.chdir('{cwd}');".format(cwd=os.getcwd())

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
    container_img,
    singularity_args,
    env_modules,
    bench_record,
    jobid,
    bench_iteration,
    cleanup_scripts,
    shadow_dir,
    edit=None,
):
    """
    Load a script from the given basedir + path and execute it.
    """
    draft = False
    if edit is not None:
        if urlparse(path).scheme == "":
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
        path, source, language = get_source(path, basedir)
    else:
        source = None

    exec_class = get_exec_class(language)

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

    if draft:
        executor.draft(listen=edit)
    else:
        executor.evaluate(edit=edit)
