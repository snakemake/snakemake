import subprocess
import shutil
import os
from urllib.parse import urlparse
import hashlib
from distutils.version import LooseVersion

import snakemake
from snakemake import conda
from snakemake.common import lazy_property, SNAKEMAKE_SEARCHPATH
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger


SNAKEMAKE_MOUNTPOINT = "/mnt/snakemake"


class Image:
    def __init__(self, url, dag):
        if " " in url:
            raise WorkflowError("Invalid singularity image URL containing "
                                "whitespace.")

        if not shutil.which("singularity"):
            raise WorkflowError("The singularity command has to be "
                                "available in order to use singularity "
                                "integration.")
        try:
            v = subprocess.check_output(["singularity", "--version"],
                                        stderr=subprocess.PIPE).decode()
        except subprocess.CalledProcessError as e:
            raise WorkflowError(
                "Failed to get singularity version:\n{}".format(
                    e.stderr.decode()))
        v = v.rsplit(" ", 1)[-1]
        if not LooseVersion(v) >= LooseVersion("2.4.1"):
            raise WorkflowError("Minimum singularity version is 2.4.1.")

        self.url = url
        self._img_dir = dag.workflow.persistence.singularity_img_path

    @property
    def is_local(self):
        scheme = urlparse(self.url).scheme
        return not scheme or scheme == "file"

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
                p = subprocess.check_output(["singularity", "pull",
                    "--name", "{}.simg".format(self.hash), self.url],
                    cwd=self._img_dir,
                    stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                raise WorkflowError("Failed to pull singularity image "
                                    "from {}:\n{}".format(self.url,
                                                          e.stdout.decode()))

    @property
    def path(self):
        if self.is_local:
            return urlparse(self.url).path
        return os.path.join(self._img_dir, self.hash) + ".simg"


def shellcmd(img_path, cmd, args="", envvars=None):
    """Execute shell command inside singularity container given optional args
       and environment variables to be passed."""

    if envvars:
        envvars = " ".join("SINGULARITYENV_{}={}".format(k, v)
                           for k, v in envvars.items())
    else:
        envvars = ""

    # mount host snakemake module into container
    args += " --bind {}:{}".format(SNAKEMAKE_SEARCHPATH, SNAKEMAKE_MOUNTPOINT)

    cmd = "{} singularity exec --home {} {} {} bash -c '{}'".format(
        envvars, os.getcwd(), args, img_path, cmd.replace("'", r"'\''"))
    return cmd
