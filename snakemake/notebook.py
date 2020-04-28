import os, sys
from urllib.error import URLError
import tempfile
import re

from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.script import get_source, ScriptBase, PythonScript, RScript

KERNEL_STARTED_RE = re.compile("Kernel started: (?P<kernel_id>\S+)")
KERNEL_SHUTDOWN_RE = re.compile("Kernel shutdown: (?P<kernel_id>\S+)")


class JupyterNotebook(ScriptBase):

    editable = True

    def draft(self):
        preamble = self.get_preamble()
        nb = nbformat.v4.new_notebook()
        self.insert_preamble_cell(preamble, nb)

        with open(self.local_path, "wb") as out:
            out.write(nbformat.writes(nb).encode())

        self.evaluate(edit=True)

    def write_script(self, preamble, fd):
        import nbformat

        nb = nbformat.reads(self.source, as_version=4)  # nbformat.NO_CONVERT
        self.insert_preamble_cell(preamble, nb)

        fd.write(nbformat.writes(nb).encode())

    def execute_script(self, fname, edit=False):
        fname_out = self.log.get("notebook", None)
        if fname_out is None or edit:
            output_parameter = ""
        else:
            fname_out = os.path.join(os.getcwd(), fname_out)
            output_parameter = "--output {fname_out:q}"

        if edit:
            cmd = "jupyter notebook --NotebookApp.quit_button=True {{fname:q}}"
        else:
            cmd = "jupyter-nbconvert --execute {output_parameter} --to notebook --ExecutePreprocessor.timeout=-1 {{fname:q}}".format(
                output_parameter=output_parameter
            )

        self._execute_cmd(cmd, fname_out=fname_out, fname=fname)

        if edit:
            if open(fname).read() == open(self.local_path).read():
                logger.warning(
                    "Notebook was not changed upon editing. Did you forget to save?"
                )
            # copy saved content back into path
            shutil.copyfile(fname, self.local_path)

    def insert_preamble_cell(self, preamble, notebook):
        preamble_cell = nbformat.v4.new_code_cell(preamble)
        preamble_cell.metadata.tags = ["snakemake-job-properties"]
        notebook["cells"].insert(0, preamble_cell)


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


def import_notebook_from_stdin(path):
    import nbformat

    nb = nbformat.read(sys.stdin, as_version=4)
    preambles = [
        i
        for i, cell in enumerate(nb["cells"])
        if "snakemake-job-properties" in cell.metadata.tags
    ]
    if len(preambles) > 1:
        raise WorkflowError(
            "More than one snakemake preamble cell found in notebook. "
            "Please clean up the notebook first, by removing all or all but one of them."
        )
    elif len(preambles) == 1:
        preamble = preambles[0]
        # remove old preamble
        del nb["cells"][preamble]

    with open(path, "wb") as out:
        out.write(nbformat.writes(nb).encode())


def get_exec_class(language):
    ExecClass = {
        "jupyter_python": PythonJupyterNotebook,
        "jupyter_r": RJupyterNotebook,
    }.get(language, None)
    if ExecClass is None:
        raise ValueError("Unsupported notebook: Expecting Jupyter Notebook (.ipynb).")
    return ExecClass


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
    edit=False,
):
    """
    Load a script from the given basedir + path and execute it.
    """
    draft = False
    if edit:
        if urlparse(path).scheme != "":
            if not os.path.exists(path):
                # draft the notebook, it does not exist yet
                language = None
                draft = True
                if path.endswith(".py.ipynb"):
                    language = "python"
                elif path.endswith(".r.ipynb"):
                    language = "r"
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

    ExecClass = get_exec_class(language)

    executor = ExecClass(
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
        executor.draft()
    else:
        executor.evaluate(edit=edit)
