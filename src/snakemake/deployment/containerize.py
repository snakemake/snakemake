from pathlib import Path
import hashlib
import os
from abc import ABC, abstractmethod

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.sourcecache import LocalSourceFile, infer_source_file
from snakemake.io import is_callable, contains_wildcard

CONDA_ENV_PATH = "/conda-envs"


# ── Step 1: The ABC (Abstract Base Class) ──
# This is the "form" that every container format must fill in.
# It defines WHAT information to output, but not HOW.


class ContainerFormat(ABC):
    """Base class for container definition formats."""

    @abstractmethod
    def header(self, image, env_hash):
        """Print the base image and metadata header."""
        pass

    @abstractmethod
    def label(self, key, value):
        """Print a metadata label."""
        pass

    @abstractmethod
    def comment(self, text):
        """Print a comment line."""
        pass

    @abstractmethod
    def copy_file(self, src, dest):
        """Print a copy/add file instruction."""
        pass

    @abstractmethod
    def add_remote_file(self, url, dest):
        """Print an instruction to download a remote file."""
        pass

    @abstractmethod
    def run_command(self, cmd):
        """Print a mkdir or similar single command."""
        pass

    @abstractmethod
    def run_combined(self, commands, final_cmd):
        """Print multiple commands combined together (e.g., conda env creation)."""
        pass


# ── Step 2: Docker format ──
# Fills in the "form" with Dockerfile syntax.
# This produces the SAME output as the original code.


class DockerFormat(ContainerFormat):

    def header(self, image, env_hash):
        print(f"FROM {image}")
        print('LABEL io.github.snakemake.containerized="true"')
        print(f'LABEL io.github.snakemake.conda_env_hash="{env_hash}"')

    def label(self, key, value):
        print(f'LABEL {key}="{value}"')

    def comment(self, text):
        print(f"# {text}")

    def copy_file(self, src, dest):
        print(f"COPY {src} {dest}")

    def add_remote_file(self, url, dest):
        print(f"ADD {url} {dest}")

    def run_command(self, cmd):
        print(f"RUN {cmd}")

    def run_combined(self, commands, final_cmd):
        all_cmds = commands + [final_cmd]
        print("\nRUN", " \\\n    ".join(all_cmds))


# ── Step 3: Apptainer format ──
# Fills in the same "form" but with Apptainer definition syntax.


class ApptainerFormat(ContainerFormat):

    def __init__(self):
        # Apptainer groups things into sections, so we collect them first
        self._labels = []
        self._files = []
        self._post = []

    def header(self, image, env_hash):
        print("Bootstrap: docker")
        print(f"From: {image}")
        self._labels.append(("io.github.snakemake.containerized", "true"))
        self._labels.append(("io.github.snakemake.conda_env_hash", env_hash))

    def label(self, key, value):
        self._labels.append((key, value))

    def comment(self, text):
        # Comments in Apptainer use #, same as Docker
        print(f"# {text}")

    def copy_file(self, src, dest):
        self._files.append((src, dest))

    def add_remote_file(self, url, dest):
        # Apptainer doesn't have ADD, so we download in %post
        self._post.append(f"mkdir -p $(dirname {dest})")
        self._post.append(f"curl -L '{url}' > '{dest}'")

    def run_command(self, cmd):
        self._post.append(cmd)

    def run_combined(self, commands, final_cmd):
        for cmd in commands:
            # Remove trailing " &&" since Apptainer runs each line separately
            self._post.append(cmd.rstrip(" &"))
        self._post.append(final_cmd)

    def finish(self):
        """Print all collected sections. Called at the end."""
        if self._labels:
            print("\n%labels")
            for key, value in self._labels:
                print(f"    {key} {value}")

        if self._files:
            print("\n%files")
            for src, dest in self._files:
                print(f"    {src} {dest}")

        if self._post:
            print("\n%post")
            for cmd in self._post:
                print(f"    {cmd}")


# ── Step 4: The main function (refactored) ──
# Lines 13–72 from the original are UNCHANGED.
# Lines 74–112 now use the format class instead of hardcoded print().


def containerize(workflow, dag, fmt="dockerfile"):
    """
    Creates a container definition for the workflow.

    Arguments:
        workflow: The workflow object.
        dag: The DAG object.
        fmt: The output format — "dockerfile" or "apptainer".
    """
    from .conda import CondaEnvFileSpec, CondaEnvSpecType

    # ── Choose the format ──
    if fmt == "dockerfile":
        formatter = DockerFormat()
    elif fmt == "apptainer":
        formatter = ApptainerFormat()
    else:
        raise WorkflowError(f"Unknown container format: {fmt}")

    # ── Everything below is the SAME as the original (lines 22–72) ──

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
        envhash.update(env.content)

    # ── From here, we use the formatter instead of hardcoded print() ──

    # Step 1: Header (was lines 74-76)
    formatter.header("condaforge/miniforge3:latest", envhash.hexdigest())

    # Step 2: Retrieve conda environments (was lines 78-109)
    generate_env_cmds = []
    generated = set()

    # curl is needed for Apptainer whenever remote env files are fetched in %post
    needs_curl = fmt == "apptainer" and any(
        not isinstance(env.file, LocalSourceFile) for env in envs
    )
    if needs_curl:
        formatter.run_command("apt-get update")
        formatter.run_command("apt-get install -y --no-install-recommends curl")

    for env in envs:
        if env.content_hash in generated:
            continue
        prefix = Path(CONDA_ENV_PATH) / env.content_hash
        env_source_path = relfile(env)
        env_target_path = prefix / "environment.yaml"

        formatter.comment(f"Conda environment:")
        formatter.comment(f"  source: {env_source_path}")
        formatter.comment(f"  prefix: {prefix}")
        for line in env.content.decode().strip().split("\n"):
            formatter.comment(f"  {line}")

        formatter.run_command(f"mkdir -p {prefix}")

        if isinstance(env.file, LocalSourceFile):
            formatter.copy_file(env_source_path, env_target_path)
        else:
            formatter.add_remote_file(
                env.file.get_path_or_uri(secret_free=True), env_target_path
            )

        generate_env_cmds.append(
            f"conda env create --prefix {prefix} --file {env_target_path} &&"
        )
        generated.add(env.content_hash)

    # Step 3: Generate conda environments (was lines 111-112)
    if generate_env_cmds:
        formatter.run_combined(generate_env_cmds, "conda clean --all -y")

    # Remove lists to reduce image size
    if needs_curl:
        formatter.run_command("rm -rf /var/lib/apt/lists/*")

    # For Apptainer, print the collected sections at the end
    if hasattr(formatter, "finish"):
        formatter.finish()
