__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import subprocess
import shutil
from distutils.version import StrictVersion
import json
from glob import glob
import tarfile
import zipfile
import uuid
from enum import Enum
import threading

from snakemake.exceptions import CreateEnvironmentException, WorkflowError
from snakemake.logging import logger
from snakemake.common import parse_uri, ON_WINDOWS
from snakemake.deployment import singularity, containerize
from .base import EnvBase, EnvCommandBase


class CondaCleanupMode(Enum):
    tarballs = "tarballs"
    cache = "cache"

    def __str__(self):
        return self.value


class Env(EnvBase):

    """Conda environment from a given specification file."""

    def __init__(
        self, env_file, workflow, env_dir=None, container_img=None, cleanup=None
    ):
        self.name = "conda"
        super().__init__(env_file, workflow, cleanup)
        self.frontend = workflow.conda_frontend
        self._container_img = container_img
        self._env_dir = env_dir or (
            containerize.CONDA_ENV_PATH
            if self.is_containerized
            else workflow.persistence.conda_env_path
        )

    @property
    def _env_archive_dir(self):
        return self.workflow.persistence.conda_env_archive_path

    def create_archive(self):
        """Create self-contained archive of environment."""
        from snakemake.shell import shell

        try:
            import yaml
        except ImportError:
            raise WorkflowError(
                "Error importing PyYAML. " "Please install PyYAML to archive workflows."
            )
        # importing requests locally because it interferes with instantiating conda environments
        import requests

        env_archive = self.archive_file
        if os.path.exists(env_archive):
            return env_archive

        try:
            # Download
            logger.info(
                "Downloading packages for conda environment {}...".format(
                    self.file.get_path_or_uri()
                )
            )
            os.makedirs(env_archive, exist_ok=True)
            try:
                out = shell.check_output(
                    "conda list --explicit --prefix '{}'".format(self.path),
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                )
                logger.debug(out)
            except subprocess.CalledProcessError as e:
                raise WorkflowError("Error exporting conda packages:\n" + e.output)
            with open(os.path.join(env_archive, "packages.txt"), "w") as pkg_list:
                for l in out.split("\n"):
                    if l and not l.startswith("#") and not l.startswith("@"):
                        pkg_url = l
                        logger.info(pkg_url)
                        parsed = parse_uri(pkg_url)
                        pkg_name = os.path.basename(parsed.uri_path)
                        # write package name to list
                        print(pkg_name, file=pkg_list)
                        # download package
                        pkg_path = os.path.join(env_archive, pkg_name)
                        with open(pkg_path, "wb") as copy:
                            r = requests.get(pkg_url)
                            r.raise_for_status()
                            copy.write(r.content)
                        try:
                            if pkg_path.endswith(".conda"):
                                assert zipfile.ZipFile(pkg_path).testzip() is None
                            else:
                                tarfile.open(pkg_path)
                        except:
                            raise WorkflowError(
                                "Package is invalid tar/zip archive: {}".format(pkg_url)
                            )
        except (
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.HTTPError,
        ):
            shutil.rmtree(env_archive)
            raise WorkflowError("Error downloading conda package {}.".format(pkg_url))
        except (Exception, BaseException) as e:
            shutil.rmtree(env_archive)
            raise e
        return env_archive

    def _create(self, env_path, env_file):
        """Creation logic for conda environment, called by self.create"""
        from snakemake.shell import shell

        # Check if env archive exists. Use that if present.
        env_archive = self.archive_file
        if os.path.exists(env_archive):
            logger.info("Installing archived conda packages.")
            pkg_list = os.path.join(env_archive, "packages.txt")
            if os.path.exists(pkg_list):
                # read pacakges in correct order
                # this is for newer env archives where the package list
                # was stored
                packages = [
                    os.path.join(env_archive, pkg.rstrip()) for pkg in open(pkg_list)
                ]
            else:
                # guess order
                packages = glob(os.path.join(env_archive, "*.tar.bz2"))

            # install packages manually from env archive
            cmd = " ".join(
                [
                    "conda",
                    "create",
                    "--quiet",
                    "--yes",
                    "--prefix '{}'".format(env_path),
                ]
                + packages
            )
            if self._container_img:
                cmd = singularity.shellcmd(
                    self._container_img.path,
                    cmd,
                    args=self._singularity_args,
                    envvars=self.get_singularity_envvars(),
                )
            out = shell.check_output(
                cmd, stderr=subprocess.STDOUT, universal_newlines=True
            )

        else:
            # Copy env file to env_path (because they can be on
            # different volumes and singularity should only mount one).
            # In addition, this allows to immediately see what an
            # environment in .snakemake/conda contains.
            target_env_file = env_path + ".yaml"
            shutil.copy(env_file, target_env_file)

            logger.info("Downloading and installing remote packages.")
            cmd = " ".join(
                [
                    self.frontend,
                    "env",
                    "create",
                    "--quiet",
                    '--file "{}"'.format(target_env_file),
                    '--prefix "{}"'.format(env_path),
                ]
            )
            if self._container_img:
                cmd = singularity.shellcmd(
                    self._container_img.path,
                    cmd,
                    args=self._singularity_args,
                    envvars=self.get_singularity_envvars(),
                )
            out = shell.check_output(
                cmd, stderr=subprocess.STDOUT, universal_newlines=True
            )

            # cleanup if requested
            if self._cleanup is CondaCleanupMode.tarballs:
                logger.info("Cleaning up conda package tarballs.")
                shell.check_output("conda clean -y --tarballs")
            elif self._cleanup is CondaCleanupMode.cache:
                logger.info("Cleaning up conda package tarballs and package cache.")
                shell.check_output("conda clean -y --tarballs --packages")

            logger.debug(out)

    @classmethod
    def get_singularity_envvars(cls):
        return {"CONDA_PKGS_DIRS": "/tmp/conda/{}".format(uuid.uuid4())}


class Conda(EnvCommandBase):
    instances = dict()
    lock = threading.Lock()
    name = "conda"

    def __new__(cls, container_img=None, prefix_path=None):
        with cls.lock:
            if container_img not in cls.instances:
                inst = super().__new__(cls)
                cls.instances[container_img] = inst
                return inst
            else:
                return cls.instances[container_img]

    def __init__(self, container_img=None, prefix_path=None):
        if not self.is_initialized:  # avoid superfluous init calls
            from snakemake.deployment import singularity
            from snakemake.shell import shell

            if isinstance(container_img, singularity.Image):
                container_img = container_img.path
            self.container_img = container_img

            if prefix_path is None or container_img is not None:
                self.prefix_path = json.loads(
                    shell.check_output(self._get_cmd("conda info --json"))
                )["conda_prefix"]
            else:
                self.prefix_path = prefix_path

            super().__init__()

    @property
    def is_initialized(self):
        return hasattr(self, "prefix_path")

    def _get_cmd(self, cmd):
        if self.container_img:
            return singularity.shellcmd(self.container_img, cmd, quiet=True)
        return cmd

    def _check(self):
        """
        Additional checks for conda (for version, etc.)
        """
        from snakemake.shell import shell

        try:
            version = shell.check_output(
                self._get_cmd("conda --version"),
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
            if self.container_img:
                version = "\n".join(
                    filter(
                        lambda line: not line.startswith("WARNING:")
                        and not line.startswith("ERROR:"),
                        version.splitlines(),
                    )
                )

            version = version.split()[1]
            if StrictVersion(version) < StrictVersion("4.2"):
                raise CreateEnvironmentException(
                    "Conda must be version 4.2 or later, found version {}.".format(
                        version
                    )
                )
        except subprocess.CalledProcessError as e:
            raise CreateEnvironmentException(
                "Unable to check conda version:\n" + e.output.decode()
            )

    def bin_path(self):
        if ON_WINDOWS:
            return os.path.join(self.prefix_path, "Scripts")
        else:
            return os.path.join(self.prefix_path, "bin")

    def shellcmd(self, env_path, cmd):
        # get path to activate script
        activate = os.path.join(self.bin_path(), "activate")

        if ON_WINDOWS:
            activate = activate.replace("\\", "/")
            env_path = env_path.replace("\\", "/")

        return "source {} '{}'; {}".format(activate, env_path, cmd)

    def shellcmd_win(self, env_path, cmd):
        """Prepend the windows activate bat script."""
        # get path to activate script
        activate = os.path.join(self.bin_path(), "activate.bat").replace("\\", "/")
        env_path = env_path.replace("\\", "/")

        return '"{}" "{}"&&{}'.format(activate, env_path, cmd)


def is_mamba_available():
    return shutil.which("mamba") is not None
