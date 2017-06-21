import os
import subprocess
import tempfile
from urllib.request import urlopen
from urllib.parse import urlparse
import hashlib
import shutil
from distutils.version import StrictVersion
import json
from glob import glob

from snakemake.exceptions import CreateCondaEnvironmentException, WorkflowError
from snakemake.logging import logger
from snakemake.common import strip_prefix
from snakemake import utils


class Env:

    """Conda environment from a given specification file."""

    def __init__(self, env_file, dag):
        self.file = env_file

        self._env_dir = dag.workflow.persistence.conda_env_path
        self._env_archive_dir = dag.workflow.persistence.conda_env_archive_path

        self._hash = None
        self._content = None
        self._path = None
        self._archive_file = None

    @property
    def content(self):
        if self._content is None:
            env_file = self.file
            if urlparse(env_file).scheme:
                content = urlopen(env_file).read()
            else:
                with open(env_file, 'rb') as f:
                    content = f.read()
            self._content = content
        return self._content

    @property
    def hash(self):
        if self._hash is None:
            md5hash = hashlib.md5()
            # Include the absolute path of the target env dir into the hash.
            # By this, moving the working directory around automatically
            # invalidates all environments. This is necessary, because binaries
            # in conda environments can contain hardcoded absolute RPATHs.
            assert os.path.isabs(self._env_dir)
            md5hash.update(self._env_dir.encode())
            md5hash.update(self.content)
            self._hash = md5hash.hexdigest()
        return self._hash

    @property
    def path(self):
        """Path to directory of the conda environment.

        First tries full hash, if it does not exist, (8-prefix) is used
        as default.

        """
        hash = self.hash
        env_dir = self._env_dir
        for h in [hash, hash[:8]]:
            path = os.path.join(env_dir, h)
            if os.path.exists(path):
                return path
        return path

    @property
    def archive_file(self):
        """Path to archive of the conda environment, which may or may not exist."""
        if self._archive_file is None:
            self._archive_file = os.path.join(self._env_archive_dir, self.hash)
        return self._archive_file

    def create_archive(self):
        """Create self-contained archive of environment."""
        try:
            import yaml
        except ImportError:
            raise WorkflowError("Error importing PyYAML. "
                "Please install PyYAML to archive workflows.")
        # importing requests locally because it interferes with instantiating conda environments
        import requests

        env_archive = self.archive_file
        if os.path.exists(env_archive):
            return env_archive

        try:
            # Download
            logger.info("Downloading packages for conda environment {}...".format(self.file))
            os.makedirs(env_archive, exist_ok=True)
            try:
                out = subprocess.check_output(["conda", "list", "--explicit",
                    "--prefix", self.path],
                    stderr=subprocess.STDOUT)
                logger.debug(out.decode())
            except subprocess.CalledProcessError as e:
                raise WorkflowError("Error exporting conda packages:\n" +
                                    e.output.decode())
            for l in out.decode().split("\n"):
                if l and not l.startswith("#") and not l.startswith("@"):
                    pkg_url = l
                    logger.info(pkg_url)
                    parsed = urlparse(pkg_url)
                    pkg_name = os.path.basename(parsed.path)
                    with open(os.path.join(env_archive, pkg_name), "wb") as copy:
                        copy.write(requests.get(pkg_url).content)
        except (Exception, BaseException) as e:
            shutil.rmtree(env_archive)
            raise e
        return env_archive

    def create(self, dryrun=False):
        """ Create the conda enviroment."""
        # Read env file and create hash.
        env_file = self.file
        tmp_file = None

        url_scheme, *_ = urlparse(env_file)
        if url_scheme and not url_scheme == 'file':
            with tempfile.NamedTemporaryFile(delete=False, suffix=".yaml") as tmp:
                tmp.write(self.content)
                env_file = tmp.name
                tmp_file = tmp.name

        env_hash = self.hash
        env_path = self.path
        # Create environment if not already present.
        if not os.path.exists(env_path):
            if dryrun:
                logger.info("Conda environment {} will be created.".format(utils.simplify_path(self.file)))
                return env_path
            logger.info("Creating conda environment {}...".format(
                        utils.simplify_path(self.file)))
            # Check if env archive exists. Use that if present.
            env_archive = self.archive_file
            try:
                if os.path.exists(env_archive):
                    # install packages manually from env archive
                    out = subprocess.check_output(["conda", "create", "--copy", "--prefix", env_path] +
                        glob(os.path.join(env_archive, "*.tar.bz2")),
                        stderr=subprocess.STDOUT
                    )
                else:
                    out = subprocess.check_output(["conda", "env", "create",
                                                "--file", env_file,
                                                "--prefix", env_path],
                                                stderr=subprocess.STDOUT)
                logger.debug(out.decode())
                logger.info("Environment for {} created (location: {})".format(
                            os.path.relpath(env_file), os.path.relpath(env_path)))
            except subprocess.CalledProcessError as e:
                # remove potential partially installed environment
                shutil.rmtree(env_path, ignore_errors=True)
                raise CreateCondaEnvironmentException(
                    "Could not create conda environment from {}:\n".format(env_file) +
                    e.output.decode())

        if tmp_file:
            # temporary file was created
            os.remove(tmp_file)

        return env_path


def check_conda():
    if shutil.which("conda") is None:
        raise CreateCondaEnvironmentException("The 'conda' command is not available in $PATH.")
    try:
        version = subprocess.check_output(["conda", "--version"], stderr=subprocess.STDOUT).decode().split()[1]
        if StrictVersion(version) < StrictVersion("4.2"):
            raise CreateCondaEnvironmentException(
                "Conda must be version 4.2 or later."
            )
    except subprocess.CalledProcessError as e:
        raise CreateCondaEnvironmentException(
            "Unable to check conda version:\n" + e.output.decode()
        )
