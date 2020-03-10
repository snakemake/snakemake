import os

from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from snakemake.script import get_source, ScriptBase, PythonScript, RScript


class JupyterNotebook(ScriptBase):
    def write_script(self, preamble, fd):
        import nbformat

        nb = nbformat.reads(self.source, as_version=4)  # nbformat.NO_CONVERT

        preamble_cell = nbformat.v4.new_code_cell(preamble)
        nb["cells"].insert(0, preamble_cell)

        fd.write(nbformat.writes(nb).encode())

    def execute_script(self, fname):
        fname_out = self.log.get("notebook", None)
        if fname_out is None:
            output_parameter = ""
        else:
            fname_out = os.path.join(os.getcwd(), fname_out)
            output_parameter = "--output {fname_out:q}"

        cmd_tmp = "jupyter-nbconvert --execute {output_parameter} --to notebook --ExecutePreprocessor.timeout=-1 {{fname:q}}".format(
            output_parameter=output_parameter
        )

        self._execute_cmd(cmd_tmp, fname_out=fname_out, fname=fname)


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
):
    """
    Load a script from the given basedir + path and execute it.
    """
    path, source, language = get_source(path, basedir)

    ExecClass = {
        "jupyter_python": PythonJupyterNotebook,
        "jupyter_r": RJupyterNotebook,
    }.get(language, None)
    if ExecClass is None:
        raise ValueError("Unsupported notebook: Expecting Jupyter Notebook (.ipynb).")

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
    executor.evaluate()
