from snakemake.sourcecache import infer_source_file
from snakemake.common import is_local_file
from pathlib import Path
import hashlib
import os
from abc import ABC, abstractmethod

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.sourcecache import LocalSourceFile
from snakemake_software_deployment_plugin_conda import EnvSpec as CondaEnvSpec

# TODO convert to rattler or pixi?


def get_containerized_path(env: CondaEnvSpec) -> Path:
    return Path("/conda-envs") / env.hash()


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

    # ── Choose the format ──
    if fmt == "dockerfile":
        formatter = DockerFormat()
    elif fmt == "apptainer":
        formatter = ApptainerFormat()
    else:
        raise WorkflowError(f"Unknown container format: {fmt}")

    # ── Everything below is the SAME as the original (lines 22–72) ──

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

    def relfile(path_or_uri):
        if is_local_file(path_or_uri):
            return os.path.relpath(path_or_uri, os.getcwd())
        else:
            return infer_source_file(path_or_uri).get_path_or_uri(secret_free=True)

    sorted_envs = sorted(conda_envs, key=lambda env: env.spec.envfile)
    envhash = hashlib.sha256()
    for env in sorted_envs:
        logger.info(f"Hashing conda environment {env.spec}.")
        # build a hash of the environment contents
        envhash.update(env.hash())

    # ── From here, we use the formatter instead of hardcoded print() ──

    # Step 1: Header (was lines 74-76)
    formatter.header("condaforge/miniforge3:latest", envhash.hexdigest())

    # Step 2: Retrieve conda environments (was lines 78-109)
    generate_env_cmds = []
    generated = set()
    for env in sorted_envs:
        if env.hash() in generated:
            # another conda env with the same content was generated before
            continue
        prefix = get_containerized_path(env)
        env_source_path = relfile(env.spec.envfile.path_or_uri)
        env_target_path = prefix / "environment.yaml"
        with open(env_source_path, "r") as f:
            env_content = f.read()
        formatter.comment(f"Conda environment:")
        formatter.comment(f"  source: {env_source_path}")
        formatter.comment(f"  prefix: {prefix}")
        for line in env_content.strip().split("\n"):
            formatter.comment(f"  {line}")
        formatter.run_command(f"mkdir -p {prefix}")

        if is_local_file(env.spec.envfile.path_or_uri):
            formatter.copy_file(env.spec.envfile.path_or_uri, env_target_path)
        else:
            formatter.add_remote_file(env.spec.envfile.path_or_uri, env_target_path)

        generate_env_cmds.append(
            f"conda env create --prefix {prefix} --file {env_target_path} &&"
        )
        generated.add(env.hash())

    # Step 3: Generate conda environments (was lines 111-112)
    if generate_env_cmds:
        formatter.run_combined(generate_env_cmds, "conda clean --all -y")

    # For Apptainer, print the collected sections at the end
    if hasattr(formatter, "finish"):
        formatter.finish()
