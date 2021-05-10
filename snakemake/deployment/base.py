__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import re
import subprocess
import tempfile
from urllib.request import urlopen
from urllib.parse import urlparse
from urllib.error import URLError
import hashlib
import shutil
from distutils.version import StrictVersion
import json
from glob import glob
import tarfile
import zipfile
import uuid
from enum import Enum
import threading
import shutil

from snakemake.exceptions import CreateCondaEnvironmentException, WorkflowError
from snakemake.logging import logger
from snakemake.common import strip_prefix, ON_WINDOWS
from snakemake import utils
from snakemake.deployment import singularity, containerize
from snakemake.io import git_content


class CleanupMode(Enum):
    tarballs = "tarballs"
    cache = "cache"

    def __str__(self):
        return self.value


class EnvBase:

    """General environment from a given specification file."""

    @property
    def content(self):
        if self._content is None:
            self._content = self._get_content()
        return self._content

    def _get_content(self):
        return self.workflow.sourcecache.open(self.file, "rb").read()

    @property
    def content_hash(self):
        if self._content_hash is None:
            md5hash = hashlib.md5()
            md5hash.update(self.content)
            self._content_hash = md5hash.hexdigest()
        return self._content_hash

    @property
    def archive_file(self):
        """Path to archive of the conda environment, which may or may not exist."""
        if self._archive_file is None:
            self._archive_file = os.path.join(self._env_archive_dir, self.content_hash)
        return self._archive_file

    def check_broken_environment(self, env_path, dryrun=False):
        """Given an environment path, check if it's broken.

        Broken means that the env_setup_start exists, but not env_setup_done.
        We resolve by removing the environment path to start fresh.
        """
        if os.path.exists(
            os.path.join(env_path, "env_setup_start")
        ) and not os.path.exists(os.path.join(env_path, "env_setup_done")):
            if dryrun:
                logger.info(
                    "Incomplete environment {} will be recreated.".format(
                        utils.simplify_path(self.file)
                    )
                )
            else:
                logger.info(
                    "Removing incomplete environment {}...".format(
                        utils.simplify_path(self.file)
                    )
                )
                shutil.rmtree(env_path, ignore_errors=True)

    @classmethod
    def get_singularity_envvars(self):
        return {"CONDA_PKGS_DIRS": "/tmp/conda/{}".format(uuid.uuid4())}

    def __hash__(self):
        # this hash is only for object comparison, not for env paths
        return hash(self.file)

    def __eq__(self, other):
        if isinstance(other, Env):
            return self.file == other.file
        return False
