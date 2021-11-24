__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import subprocess
import shutil
import os
import hashlib
from distutils.version import LooseVersion

import snakemake
from snakemake.deployment.conda import Conda
from snakemake.common import (
    is_local_file,
    parse_uri,
    lazy_property,
    SNAKEMAKE_SEARCHPATH,
)
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger


SNAKEMAKE_MOUNTPOINT = "/mnt/snakemake"


class Image:
    def __init__(self, url, dag, is_containerized):
        if " " in url:
            raise WorkflowError(
                "Invalid singularity image URL containing " "whitespace."
            )

        if not shutil.which("singularity"):
            raise WorkflowError(
                "The singularity command has to be "
                "available in order to use singularity "
                "integration."
            )
        try:
            v = subprocess.check_output(
                ["singularity", "--version"], stderr=subprocess.PIPE
            ).decode()
        except subprocess.CalledProcessError as e:
            raise WorkflowError(
                "Failed to get singularity version:\n{}".format(e.stderr.decode())
            )
        v = v.rsplit(" ", 1)[-1]
        if v.startswith("v"):
            v = v[1:]
        if not LooseVersion(v) >= LooseVersion("2.4.1"):
            raise WorkflowError("Minimum singularity version is 2.4.1.")

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
        if self.is_local:
            return
        if dryrun:
            logger.info("Singularity image {} will be pulled.".format(self.url))
            return
        logger.debug("Singularity image location: {}".format(self.path))
        if not os.path.exists(self.path):
            logger.info("Pulling singularity image {}.".format(self.url))
            try:
                p = subprocess.check_output(
                    [
                        "singularity",
                        "pull",
                        "--name",
                        "{}.simg".format(self.hash),
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
        envvars = " ".join(
            "SINGULARITYENV_{}={}".format(k, v) for k, v in envvars.items()
        )
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
        args += " --bind {}:{}".format(SNAKEMAKE_SEARCHPATH, SNAKEMAKE_MOUNTPOINT)

    if container_workdir:
        args += " --pwd {}".format(container_workdir)

    cmd = "{} singularity {} exec --home {} {} {} {} -c '{}'".format(
        envvars,
        "--quiet --silent" if quiet else "",
        os.getcwd(),
        args,
        img_path,
        shell_executable,
        cmd.replace("'", r"'\''"),
    )
    logger.debug(cmd)
    return cmd
