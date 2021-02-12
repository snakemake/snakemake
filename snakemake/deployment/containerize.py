from pathlib import Path
import hashlib
import os

from snakemake.common import get_container_image
from snakemake.io import contains_wildcard
from snakemake.exceptions import WorkflowError
from snakemake.deployment import conda
from snakemake.logging import logger


CONDA_ENV_PATH = "/conda-envs"


def containerize(workflow):
    if any(
        contains_wildcard(rule.conda_env)
        for rule in workflow.rules
        if rule.conda_env is not None
    ):
        raise WorkflowError(
            "Containerization of conda based workflows is not allowed if any conda env definition contains a wildcard."
        )

    relfile = lambda env: os.path.relpath(env.file, os.getcwd())

    envs = sorted(
        set(
            conda.Env(rule.conda_env, workflow, env_dir=CONDA_ENV_PATH)
            for rule in workflow.rules
            if rule.conda_env is not None
        ),
        key=relfile,
    )
    envhash = hashlib.sha256()
    for env in envs:
        logger.info("Hashing conda environment {}.".format(relfile(env)))
        # build a hash of the environment contents
        envhash.update(env.content)

    print("FROM condaforge/mambaforge:latest")
    print('LABEL io.github.snakemake.containerized="true"')
    print('LABEL io.github.snakemake.conda_env_hash="{}"'.format(envhash.hexdigest()))

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
        get_env_cmds.append("#   source: {}".format(env_source_path))
        get_env_cmds.append("#   prefix: {}".format(prefix))
        get_env_cmds.append(
            "\n".join(map("#   {}".format, env.content.decode().strip().split("\n")))
        )
        get_env_cmds.append("RUN mkdir -p {}".format(prefix))
        if env.file.startswith("https://"):
            # get_env_cmds.append("RUN curl -sSL {} > {}".format(env.file, env_target_path))
            get_env_cmds.append("ADD {} {}".format(env.file, env_target_path))
        else:
            get_env_cmds.append("COPY {} {}".format(env_source_path, env_target_path))

        generate_env_cmds.append(
            "mamba env create --prefix {} --file {} &&".format(prefix, env_target_path)
        )
        generated.add(env.content_hash)

    print("\n# Step 1: Retrieve conda environments")
    for cmd in get_env_cmds:
        print(cmd)

    print("\n# Step 2: Generate conda environments")
    print("\nRUN", " \\\n    ".join(generate_env_cmds), "\\\n    mamba clean --all -y")
