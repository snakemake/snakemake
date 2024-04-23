__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
import subprocess
import shutil
import os
import hashlib

from snakemake.common import (
    get_snakemake_searchpaths,
    is_local_file,
    parse_uri,
)
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake_interface_common.utils import lazy_property


SNAKEMAKE_MOUNTPOINT = "/mnt/snakemake"


def get_snakemake_searchpath_mountpoints():
    paths = get_snakemake_searchpaths()
    base = Path("/mnt/snakemake_searchpaths")
    return [str(base / f"item_{i}") for i in range(len(paths))]


class Image:
    def __init__(self, url, dag, is_containerized):
        if " " in url:
            raise WorkflowError("Invalid singularity image URL containing whitespace.")

        self.singularity = Singularity()

        self.url = url
        self._img_dir = dag.workflow.persistence.container_img_path
        self.is_containerized = is_containerized

    @property
    def is_local(self):
        return is_local_file(self.url)

    @lazy_property
    def hash(self):
        md5hash = hashlib.md5()
        md5hash.update(self.url.encode())
        return md5hash.hexdigest()

    def pull(self, dryrun=False):
        self.singularity.check()
        if self.is_local:
            return
        if dryrun:
            logger.info(f"Singularity image {self.url} will be pulled.")
            return
        logger.debug(f"Singularity image location: {self.path}")
        if not os.path.exists(self.path):
            logger.info(f"Pulling singularity image {self.url}.")
            try:
                p = subprocess.check_output(
                    [
                        "singularity",
                        "pull",
                        "--name",
                        f"{self.hash}.simg",
                        self.url,
                    ],
                    cwd=self._img_dir,
                    stderr=subprocess.STDOUT,
                )
            except subprocess.CalledProcessError as e:
                raise WorkflowError(
                    "Failed to pull singularity image "
                    "from {}:\n{}".format(self.url, e.stdout.decode())
                )

    @property
    def path(self):
        if self.is_local:
            return parse_uri(self.url).uri_path
        return os.path.join(self._img_dir, self.hash) + ".simg"

    def __hash__(self):
        return hash(self.hash)

    def __eq__(self, other):
        return self.url == other.url


def shellcmd(
    img_path,
    cmd,
    args="",
    quiet=False,
    envvars=None,
    shell_executable=None,
    container_workdir=None,
    is_python_script=False,
):
    """Execute shell command inside singularity container given optional args
    and environment variables to be passed."""

    if envvars:
        envvars = " ".join(f"SINGULARITYENV_{k}={v}" for k, v in envvars.items())
    else:
        envvars = ""

    if shell_executable is None:
        shell_executable = "sh"
    else:
        # Ensure to just use the name of the executable, not a path,
        # because we cannot be sure where it is located in the container.
        shell_executable = os.path.split(shell_executable)[-1]

    if is_python_script:
        # mount host snakemake module into container
        args += " ".join(
            f" --bind {repr(searchpath)}:{repr(mountpoint)}"
            for searchpath, mountpoint in zip(
                get_snakemake_searchpaths(), get_snakemake_searchpath_mountpoints()
            )
        )

    if container_workdir:
        args += f" --pwd {repr(container_workdir)}"

    cmd = "{} singularity {} exec --home {} {} {} {} -c '{}'".format(
        envvars,
        "--quiet --silent" if quiet else "",
        repr(os.getcwd()),
        args,
        img_path,
        shell_executable,
        cmd.replace("'", r"'\''"),
    )
    logger.debug(cmd)
    return cmd


class Singularity:
    instance = None

    def __new__(cls):
        if cls.instance is not None:
            return cls.instance
        else:
            inst = super().__new__(cls)
            cls.instance = inst
            return inst

    def __init__(self):
        self.checked = False
        self._version = None

    @property
    def version(self):
        assert (
            self._version is not None
        ), "bug: singularity version accessed before check() has been called"
        return self._version

    def parseversion(self, raw_version):
        import packaging

        raw_version = raw_version.rsplit(" ", 1)[-1]
        if raw_version.startswith("v"):
            raw_version = raw_version[1:]

        parsed_version = None
        trimend = len(raw_version)
        while parsed_version is None:
            try:
                parsed_version = packaging.version.Version(raw_version[:trimend])
            except packaging.version.InvalidVersion:
                trimend = trimend - 1
                if trimend == 0:
                    raise WorkflowError(
                        f"Apptainer/Singularity version cannot be parsed: {raw_version}"
                    )
        return parsed_version

    def check(self):
        from packaging.version import parse

        if not self.checked:
            if not shutil.which("singularity"):
                raise WorkflowError(
                    "The apptainer or singularity command has to be "
                    "available in order to use apptainer/singularity "
                    "integration."
                )
            try:
                v = subprocess.check_output(
                    ["singularity", "--version"], stderr=subprocess.PIPE
                ).decode()
            except subprocess.CalledProcessError as e:
                raise WorkflowError(
                    f"Failed to get singularity version:\n{e.stderr.decode()}"
                )
            if v.startswith("apptainer"):
                if self.parseversion(v) < parse("1.0.0"):
                    raise WorkflowError("Minimum apptainer version is 1.0.0.")
            elif self.parseversion(v) < parse("2.4.1"):
                raise WorkflowError("Minimum singularity version is 2.4.1.")
            self._version = v
