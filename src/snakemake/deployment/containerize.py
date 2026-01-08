from pathlib import Path
import hashlib
import os

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.sourcecache import LocalSourceFile, infer_source_file
from snakemake.io import is_callable, contains_wildcard


CONDA_ENV_PATH = "/conda-envs"


def containerize(workflow, dag):
    """
    Creates a Dockerfile for containerizing the workflow.

    This function iterates through both the jobs in the initial DAG and all rules
    in the workflow to collect their Conda environment specifications.
    """
    from .conda import CondaEnvFileSpec, CondaEnvSpecType

    # collect envs from jobs from the initial DAG.
    conda_envs = {job.conda_env for job in dag.jobs if job.conda_env is not None}
    dag_jobs = {job.name for job in dag.jobs if job.conda_env is not None}

    # collect envs from rules not in the initial DAG (e.g., for rules past checkpoints)
    for rule in workflow.rules:
        if not rule.conda_env or rule.name in dag_jobs:
            continue

        env_def = rule.conda_env

        if is_callable(env_def) or contains_wildcard(str(env_def)):
            raise WorkflowError(
                "Containerization of conda based workflows is not allowed if any conda env definition is dynamic "
                "(e.g. contains a wildcard or is a function) and the corresponding rule is not part of the initial DAG. "
                f"Found in rule {rule.name}."
            )

        spec_type = CondaEnvSpecType.from_spec(rule.conda_env)
        if spec_type is not CondaEnvSpecType.FILE:
            raise WorkflowError(
                "Only file-based conda environments are supported for containerization for rules that are not part "
                "of the initial DAG (e.g. because they are downstream of a checkpoint). "
                f"Rule {rule.name} uses a conda environment by name or directory."
            )

        env_def = infer_source_file(rule.conda_env, rule.basedir)
        spec = CondaEnvFileSpec(env_def)
        conda_env = spec.get_conda_env(workflow, envs_dir=CONDA_ENV_PATH)
        if conda_env.is_externally_managed:
            logger.warning(
                "Skipping containerization of externally managed conda environment in rule %s.",
                rule.name,
            )
            continue
        conda_envs.add(conda_env)

    def relfile(env):
        if isinstance(env.file, LocalSourceFile):
            return os.path.relpath(
                env.file.get_path_or_uri(secret_free=True), os.getcwd()
            )
        else:
            return env.file.get_path_or_uri(secret_free=True)

    envs = sorted(conda_envs, key=relfile)
    envhash = hashlib.sha256()
    for env in envs:
        logger.info(f"Hashing conda environment {relfile(env)}.")
        # build a hash of the environment contents
        envhash.update(env.content)

    print("FROM condaforge/miniforge3:latest")
    print('LABEL io.github.snakemake.containerized="true"')
    print(f'LABEL io.github.snakemake.conda_env_hash="{envhash.hexdigest()}"')

    generated = set()
    get_env_cmds = []
    generate_env_cmds = []
    for env in envs:
        if env.content_hash in generated:
            # another conda env with the same content was generated before
            continue
        prefix = Path(CONDA_ENV_PATH) / env.content_hash
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
