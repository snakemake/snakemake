from pathlib import Path
import hashlib
import os

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.sourcecache import LocalSourceFile, infer_source_file
from snakemake.io import is_callable, contains_wildcard
from snakemake_software_deployment_plugin_conda import EnvSpec as CondaEnvSpec

CONDA_ENV_PATH = "/conda-envs"

# TODO convert to rattler or pixi?


def containerize(workflow, dag):
    """
    Creates a Dockerfile for containerizing the workflow.

    This function iterates through both the jobs in the initial DAG and all rules
    in the workflow to collect their Conda environment specifications.
    """
    from .conda import CondaEnvFileSpec, CondaEnvSpecType

    # collect envs from jobs from the initial DAG.
    def is_conda_job(job):
        return job.software_env_spec is not None and isinstance(
            job.software_env_spec, CondaEnvSpec
        )
    conda_envs = {job.software_env for job in dag.jobs if is_conda_job(job)}

    for env in conda_envs:
        if env.spec.envfile is None:
            raise WorkflowError(
                "Only file-based conda environments are supported for containerization for rules that are not part "
                "of the initial DAG (e.g. because they are downstream of a checkpoint). "
                f"Invalid conda env spec: {env.spec}"
            )

    def relfile(env):
        if isinstance(env.file, LocalSourceFile):
            return os.path.relpath(
                env.file.get_path_or_uri(secret_free=True), os.getcwd()
            )
        else:
            return env.file.get_path_or_uri(secret_free=True)

    sorted_envs = sorted(conda_envs, key=lambda env: env.spec.envfile)
    envhash = hashlib.sha256()
    for env in sorted_envs:
        logger.info(f"Hashing conda environment {env.spec}.")
        # build a hash of the environment contents
        envhash.update(env.hash())

    print("FROM condaforge/miniforge3:latest")
    print('LABEL io.github.snakemake.containerized="true"')
    print(f'LABEL io.github.snakemake.conda_env_hash="{envhash.hexdigest()}"')

    generated = set()
    get_env_cmds = []
    generate_env_cmds = []
    for env in sorted_envs:
        if env.hash() in generated:
            # another conda env with the same content was generated before
            continue
        prefix = Path(CONDA_ENV_PATH) / env.hash()
        env_source_path = relfile(env)
        env_target_path = prefix / "environment.yaml"
        get_env_cmds.append("\n# Conda environment:")
        get_env_cmds.append(f"#   source: {env_source_path}")
        get_env_cmds.append(f"#   prefix: {prefix}")
        get_env_cmds.append(
            "\n".join(map("#   {}".format, env.content.decode().strip().split("\n")))
        )
        get_env_cmds.append(f"RUN mkdir -p {prefix}")
        if isinstance(env.file, LocalSourceFile):
            get_env_cmds.append(f"COPY {env_source_path} {env_target_path}")
        else:
            get_env_cmds.append(
                f"ADD {env.file.get_path_or_uri(secret_free=True)} {env_target_path}"
            )

        generate_env_cmds.append(
            f"conda env create --prefix {prefix} --file {env_target_path} &&"
        )
        generated.add(env.content_hash)

    print("\n# Step 2: Retrieve conda environments")
    for cmd in get_env_cmds:
        print(cmd)

    print("\n# Step 3: Generate conda environments")
    print("\nRUN", " \\\n    ".join(generate_env_cmds), "\\\n    conda clean --all -y")
