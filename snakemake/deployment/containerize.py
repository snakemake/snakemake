from pathlib import Path
import hashlib
import os

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.sourcecache import LocalSourceFile


CONDA_ENV_PATH = "/conda-envs"


def containerize(workflow, dag):
    if any(
        job.conda_env_spec.contains_wildcard
        for job in dag.jobs
        if job.conda_env_spec is not None
    ):
        raise WorkflowError(
            "Containerization of conda based workflows is not allowed if any conda env definition contains a wildcard."
        )

    def relfile(env):
        if isinstance(env.file, LocalSourceFile):
            return os.path.relpath(env.file.get_path_or_uri(), os.getcwd())
        else:
            return env.file.get_path_or_uri()

    envs = sorted(
        set(
            job.conda_env_spec.get_conda_env(workflow, envs_dir=CONDA_ENV_PATH)
            for job in dag.jobs
            if job.conda_env_spec is not None
        ),
        key=relfile,
    )
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
            get_env_cmds.append(f"ADD {env.file.get_path_or_uri()} {env_target_path}")

        generate_env_cmds.append(
            f"conda env create --prefix {prefix} --file {env_target_path} &&"
        )
        generated.add(env.content_hash)

    print("\n# Step 2: Retrieve conda environments")
    for cmd in get_env_cmds:
        print(cmd)

    print("\n# Step 3: Generate conda environments")
    print("\nRUN", " \\\n    ".join(generate_env_cmds), "\\\n    conda clean --all -y")
