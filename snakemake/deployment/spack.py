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
from .base import EnvBase

from snakemake.exceptions import CreateSpackEnvironmentException, WorkflowError
from snakemake.logging import logger
from snakemake.common import strip_prefix, ON_WINDOWS
from snakemake import utils
from snakemake.deployment import singularity, containerize
from snakemake.io import git_content


class SpackCleanupMode(Enum):
    tarballs = "tarballs"
    cache = "cache"

    def __str__(self):
        return self.value


class Env(EnvBase):

    """Spack environment from a given specification file."""

    def __init__(self, env_file, workflow, env_dir=None, cleanup=None):
        self.file = env_file
        self.workflow = workflow

        self._env_dir = env_dir or workflow.persistence.spack_env_path
        self._hash = None
        self._content_hash = None
        self._content = None
        self._path = None
        self._archive_file = None
        self._cleanup = cleanup

    @property
    def _env_archive_dir(self):
        return self.workflow.persistence.spack_env_archive_path

    @property
    def hash(self):
        if self._hash is None:
            md5hash = hashlib.md5()
            # Include the absolute path of the target env dir into the hash.
            # By this, moving the working directory around automatically
            # invalidates all environments. This is necessary, because binaries
            # in spack environments can contain hardcoded absolute RPATHs.
            env_dir = os.path.realpath(self._env_dir)
            md5hash.update(env_dir.encode())
            md5hash.update(self.content)
            self._hash = md5hash.hexdigest()
        return self._hash

    @property
    def path(self):
        """Path to directory of the spack environment.

        First tries full hash, if it does not exist, (8-prefix) is used
        as default.
        """
        hash = self.hash
        env_dir = self._env_dir
        get_path = lambda h: os.path.join(env_dir, h)
        hash_candidates = [hash[:8], hash]  # [0] is the old fallback hash (shortened)
        exists = [os.path.exists(get_path(h)) for h in hash_candidates]
        if exists[1] or (not exists[0]):
            # full hash exists or fallback hash does not exist: use full hash
            return get_path(hash_candidates[1])
        # use fallback hash
        return get_path(hash_candidates[0])

    def create_archive(self):
        """
        Create self-contained archive of environment.

        This won't be possible with spack until an environment can be created
        that uses a copy instead of link. Currently, we just include all
        of the links, the idea being that the environment could be restored
        on the same system.
        """
        from snakemake.shell import shell

        try:
            import yaml
        except ImportError:
            raise WorkflowError(
                "Error importing PyYAML. " "Please install PyYAML to archive workflows."
            )
        # importing requests locally because it interferes with instantiating spack environments
        import requests

        env_archive = self.archive_file
        if os.path.exists(env_archive):
            return env_archive

        try:
            logger.info("Saving metadata spack environment {}...".format(self.file))
            os.makedirs(env_archive, exist_ok=True)

            # Activate the environment so we can find packages
            activate = "eval `spack env activate --sh %s`" % self.path
            shell.check_output(
                activate, stderr=subprocess.STDOUT, universal_newlines=True
            )

            try:
                out = shell.check_output(
                    "spack find --json",
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                )
                logger.debug(out)
            except subprocess.CalledProcessError as e:
                raise WorkflowError("Error exporting spack spec:\n" + e.output)

            # TODO: spack has support for copy/relocate environments, but it is
            # usually the case that the paths are too long to change for the
            # relocate. So at best we could only save hard coded paths (and don't))
            with open(os.path.join(env_archive, "specs.json"), "w") as pkg_list:
                print(out, file=pkg_list)

            deactivate = "eval `spack env deactivate --sh`"
            out = shell.check_output(
                deactivate, stderr=subprocess.STDOUT, universal_newlines=True
            )

        except (Exception, BaseException) as e:
            shutil.rmtree(env_archive)
            raise e

        return env_archive

    def create(self, dryrun=False):
        """
        Create the spack enviroment.
        """
        from snakemake.shell import shell

        # Read env file and create hash.
        env_file = self.file
        tmp_file = None

        url_scheme, *_ = urlparse(env_file)
        if (url_scheme and not url_scheme == "file") or (
            not url_scheme and env_file.startswith("git+file:/")
        ):
            with tempfile.NamedTemporaryFile(delete=False, suffix=".yaml") as tmp:
                tmp.write(self.content)
                env_file = tmp.name
                tmp_file = tmp.name

        env_hash = self.hash
        env_path = self.path

        # Check for broken environment
        self.check_broken_environment(env_path, dryrun)

        # Create environment if not already present.
        if not os.path.exists(env_path):
            if dryrun:
                logger.info(
                    "Spack environment {} will be created.".format(
                        utils.simplify_path(self.file)
                    )
                )
                return env_path

            spack = Spack()
            logger.info(
                "Creating spack environment {}...".format(
                    utils.simplify_path(self.file)
                )
            )
            env_archive = self.archive_file

            # This currently doesn't work because the paths are too long to relocate
            # We need to specify it will be a copy, not a symlink
            #    view:
            #      copy:
            #        root: /Users/spackuser/soft/test/copy
            #        link_type: copy

            # spack_env_root = os.path.join(env_path, "snakemake_view")

            # with open(env_file, 'r') as fd:
            #    updated_file = yaml.load(fd.read(), Loader=yaml.loader.BaseLoader)
            # updated_file['spack']['view'] = {"copy": {"root": spack_env_root, "link_type": "copy"}}

            try:
                # Touch "start" flag file
                os.makedirs(env_path, exist_ok=True)
                with open(os.path.join(env_path, "env_setup_start"), "a") as f:
                    pass

                # TODO: check if env archive exists. Use that if present.
                # Copy env file to env_path so we can see what an
                # environment in .snakemake/spack contains.
                target_env_file = os.path.join(env_path, "spack.yaml")

                # See lines above about relocation paths
                # with open(target_env_file, 'w') as fd:
                #    yaml.dump(updated_file, fd)
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

                # Touch "done" flag file
                with open(os.path.join(env_path, "env_setup_done"), "a") as f:
                    pass

                logger.debug(out)
                logger.info(
                    "Environment for {} created (location: {})".format(
                        os.path.relpath(env_file), os.path.relpath(env_path)
                    )
                )
            except subprocess.CalledProcessError as e:
                # remove potential partially installed environment
                shutil.rmtree(env_path, ignore_errors=True)
                raise CreateSpackEnvironmentException(
                    "Could not create spack environment from {}:\n".format(env_file)
                    + e.output
                )

        if tmp_file:
            # temporary file was created
            os.remove(tmp_file)

        return env_path

    def __hash__(self):
        # this hash is only for object comparison, not for env paths
        return hash(self.file)

    def __eq__(self, other):
        if isinstance(other, Env):
            return self.file == other.file
        return False


class Spack:
    def __init__(self):
        self._check()

    def _check(self):
        from snakemake.shell import shell

        # Use type here since conda now is a function.
        # type allows to check for both functions and regular commands.
        if not ON_WINDOWS or shell.get_executable():
            locate_cmd = "type spack"
        else:
            locate_cmd = "where spack"

        try:
            shell.check_output(locate_cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            raise CreateSpackEnvironmentException(
                "The 'spack' command is not "
                "available in the "
                "shell {} that will be "
                "used by Snakemake. You have "
                "to ensure that it is in your "
                "PATH.".format(shell.get_executable())
            )

    def shellcmd(self, env_path, cmd):
        activate = "eval `spack env activate --sh %s`" % env_path
        return "{}; {}".format(activate, cmd)
