__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import subprocess
import shutil
from enum import Enum
from .base import EnvBase, EnvCommandBase

from snakemake.logging import logger


class SpackCleanupMode(Enum):
    tarballs = "tarballs"
    cache = "cache"

    def __str__(self):
        return self.value


class Env(EnvBase):

    """Spack environment from a given specification file."""

    def __init__(self, env_file, workflow, env_dir=None, cleanup=None):
        self.name = "spack"
        super().__init__(env_file, workflow, cleanup)
        self._env_dir = env_dir or workflow.persistence.spack_env_path

    @property
    def _env_archive_dir(self):
        return self.workflow.persistence.spack_env_archive_path

    def _create(self, env_path, env_file):
        """
        Create the spack environment!
        """
        from snakemake.shell import shell

        # Copy env file to env_path so we can see what an
        # environment in .snakemake/spack contains.
        target_env_file = os.path.join(env_path, "spack.yaml")
        shutil.copyfile(env_file, target_env_file)

        logger.info("Downloading and installing remote packages.")
        create = " ".join(
            [
                "spack",
                "env",
                "create",
                "--dir '{}'".format(env_path),
                target_env_file,
            ]
        )

        ## Finally, source and install (also concretized), and deactivate
        install = "spack install"
        activate = "eval `spack env activate --sh %s`; %s" % (env_path, install)

        for cmd in [create, activate]:
            logger.info(cmd)
            out = shell.check_output(
                cmd, stderr=subprocess.STDOUT, universal_newlines=True
            )
            logger.debug(out)


class Spack(EnvCommandBase):
    # name of the executable (required)
    name = "spack"

    def shellcmd(self, env_path, cmd):
        activate = "eval `spack env activate --sh %s`" % env_path
        return "{}; {}".format(activate, cmd)
