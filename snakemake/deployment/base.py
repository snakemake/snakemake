__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import snakemake
import subprocess
import tempfile
import hashlib
from enum import Enum
from pathlib import Path
import shutil

from snakemake.exceptions import CreateEnvironmentException, WorkflowError
from snakemake.logging import logger
from snakemake.deployment import singularity
from snakemake.sourcecache import LocalGitFile, LocalSourceFile, infer_source_file
from snakemake.common import ON_WINDOWS


class CleanupMode(Enum):
    tarballs = "tarballs"
    cache = "cache"

    def __str__(self):
        return self.value


class EnvCommandBase:
    """
    An environment command controller
    """

    def __init__(self):
        self.check()

    def shellcmd(self, env_path, cmd):
        raise NotImplementedError

    def check(self):
        """Every environment command should have a check to ensure requirements"""
        from snakemake.shell import shell

        if not hasattr(self, "name"):
            raise WorkflowError("An environment command is required to have a name.")

        # Use type here since conda now is a function.
        # type allows to check for both functions and regular commands.
        if not ON_WINDOWS or shell.get_executable():
            locate_cmd = "type %s" % self.name
        else:
            locate_cmd = "where %s" % self.name

        try:
            shell.check_output(locate_cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            raise CreateEnvironmentException(
                "The '{}' command is not "
                "available in the "
                "shell {} that will be "
                "used by Snakemake. You have "
                "to ensure that it is in your "
                "PATH.".format(self.name, shell.get_executable())
            )

        # Run command specific checks (or none)
        self._check()

    def _check(self):
        pass


class EnvBase:
    """
    General environment from a given specification file.
    """

    def __init__(self, env_file, workflow, cleanup=None):
        self.file = infer_source_file(env_file)
        self._cleanup = cleanup

        # Not all environments are required to support container images
        self._container_img = None
        self._singularity_args = workflow.singularity_args

        self.workflow = workflow
        self._env_dir = None
        self._hash = None
        self._content_hash = None
        self._content = None
        self._path = None
        self._archive_file = None
        self._cleanup = cleanup

        # A name is required
        if not hasattr(self, "name"):
            raise CreateEnvironmentException(
                "Environment class is missing name attribute, please add it."
            )

    @property
    def content(self):
        if self._content is None:
            self._content = self._get_content()
        return self._content

    def _get_content(self):
        return self.workflow.sourcecache.open(self.file, "rb").read()

    @property
    def container_img_url(self):
        return self._container_img.url if self._container_img else None

    @property
    def content_hash(self):
        if self._content_hash is None:
            md5hash = hashlib.md5()
            md5hash.update(self.content)
            self._content_hash = md5hash.hexdigest()
        return self._content_hash

    @property
    def hash(self):
        if self._hash is None:
            if self.is_containerized:
                self._hash = self.content_hash
            else:
                self.md5_hash()
        return self._hash

    @property
    def is_containerized(self):
        if not self._container_img:
            return False
        return self._container_img.is_containerized

    def create_archive(self):
        """
        Create self-contained archive of environment.

        This won't be possible with spack until an environment can be created
        that uses a copy instead of link (now possible) but ALSO can be
        relocated. Currently, the spack.yaml files are the best reproducers.
        """
        raise WorkflowError("%s does not support environment archive." % self.name)

    @property
    def path(self):
        """Path to directory of the environment.

        First tries full hash, if it does not exist, (8-prefix) is used
        as default.

        """
        hash = self.hash
        env_dir = self._env_dir
        get_path = lambda h: os.path.join(env_dir, h)
        hash_candidates = [hash[:8], hash]  # [0] is the old fallback hash (shortened)
        exists = [os.path.exists(get_path(h)) for h in hash_candidates]
        if self.is_containerized or exists[1] or (not exists[0]):
            # containerizes, full hash exists or fallback hash does not exist: use full hash
            return get_path(hash_candidates[1])
        # use fallback hash
        return get_path(hash_candidates[0])

    @property
    def archive_file(self):
        """Path to archive of the conda environment, which may or may not exist."""
        if self._archive_file is None:
            self._archive_file = os.path.join(self._env_archive_dir, self.content_hash)
        return self._archive_file

    def _init_env_file(self):
        """
        Create the temporary clone of the environment file and return a lookup
        """
        # Read env file and create hash.
        env_file = self.file
        tmp_file = None
        env_path = self.path

        if not isinstance(env_file, LocalSourceFile) or isinstance(
            env_file, LocalGitFile
        ):
            with tempfile.NamedTemporaryFile(delete=False, suffix=".yaml") as tmp:
                tmp.write(self.content)
                env_file = tmp.name
                tmp_file = tmp.name
        else:
            env_file = env_file.get_path_or_uri()
        return {"env_file": env_file, "tmp_file": tmp_file, "env_path": env_path}

    def md5_hash(self):
        """
        Set the self.hash to an md5hash
        """
        md5hash = hashlib.md5()
        # Include the absolute path of the target env dir into the hash.
        # By this, moving the working directory around automatically
        # invalidates all environments. This is necessary, because binaries
        # in environments can contain hardcoded absolute RPATHs.
        env_dir = os.path.realpath(self._env_dir)
        md5hash.update(env_dir.encode())
        if self._container_img:
            md5hash.update(self._container_img.url.encode())
        md5hash.update(self.content)
        self._hash = md5hash.hexdigest()

    def _create(self, env_path, env_file):
        raise NotImplementedError

    def _touch_start(self, env_path):
        self._touch_file(os.path.join(env_path, "env_setup_start"))

    def _touch_done(self, env_path):
        self._touch_file(os.path.join(env_path, "env_setup_done"))

    def _touch_file(self, touch_path):
        """
        Touch a start or done file for the environment creation process!
        """
        dirname = os.path.dirname(touch_path)
        os.makedirs(dirname, exist_ok=True)
        Path(touch_path).touch()

    def create(self, dryrun=False):
        """
        Shared create for an environment - requires a subclass _create function
        """
        from snakemake.shell import shell

        # Read env file and create hash.
        paths = self._init_env_file()
        env_file = paths["env_file"]
        tmp_file = paths["tmp_file"]
        env_path = paths["env_path"]

        # Not all environments need to support containerization, conda does
        if self.is_containerized:
            if not dryrun:
                try:
                    shell.check_output(
                        singularity.shellcmd(
                            self._container_img.path,
                            "[ -d '{}' ]".format(env_path),
                            args=self._singularity_args,
                            envvars=self.get_singularity_envvars(),
                            quiet=True,
                        ),
                        stderr=subprocess.PIPE,
                    )
                except subprocess.CalledProcessError as e:
                    raise WorkflowError(
                        "Unable to find environment in container image. "
                        "Maybe a {} environment was modified without containerizing again "
                        "(see snakemake --containerize)?\nDetails:\n{}\n{}".format(
                            self.name, e, e.stderr.decode()
                        )
                    )
                return env_path
            else:
                # env should be present in the container
                return env_path

        # Check for broken environment
        self.check_broken_environment(env_path, dryrun)

        # Create environment if not already present.
        if not os.path.exists(env_path):
            if dryrun:
                logger.info(
                    "Conda environment {} will be created.".format(
                        self.file.simplify_path()
                    )
                )
                return env_path

            logger.info(
                "Creating {} environment {}...".format(
                    self.name, self.file.simplify_path()
                )
            )
            try:
                self._touch_start(env_path)
                self._create(env_path, env_file)
                self._touch_done(env_path)
                logger.info(
                    "Environment for {} created (location: {})".format(
                        os.path.relpath(env_file), os.path.relpath(env_path)
                    )
                )

            except subprocess.CalledProcessError as e:
                # remove potential partially installed environment
                shutil.rmtree(env_path, ignore_errors=True)
                raise CreateEnvironmentException(
                    "Could not create environment from {}:\n".format(env_file)
                    + e.output
                )

        if tmp_file:
            os.remove(tmp_file)

        return env_path

    @classmethod
    def get_singularity_envvars(cls):
        return {}

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
                    "Incomplete {} environment {} will be recreated.".format(
                        self.name.capitalize(), self.file.simplify_path()
                    )
                )
            else:
                logger.info(
                    "Removing incomplete {} environment {}...".format(
                        self.name.capitalize(), self.file.simplify_path()
                    )
                )
                shutil.rmtree(env_path, ignore_errors=True)

    def __hash__(self):
        # this hash is only for object comparison, not for env paths
        return hash(self.file)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.file == other.file
        return False
