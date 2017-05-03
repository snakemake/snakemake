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


class Env:

    """Conda environment from a given specification file."""

    def __init__(self, env_file):
        self.file = env_file

    @property
    def content(self):
        if not hasattr(self, "_content"):
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
        if not hasattr(self, "_hash"):
            md5hash = hashlib.md5()
            md5hash.update(self.content)
            self._hash = md5hash.hexdigest()
        return self._hash

    def get_archive(self, dag):
        """Get path to archived environment derived from given environment file."""
        return os.path.join(dag.workflow.persistence.conda_env_archive_path, self.hash)

    def create_archive(self, dag):
        """Create self-contained archive of environment."""
        try:
            import yaml
        except ImportError:
            raise WorkflowError("Error importing PyYAML. "
                "Please install PyYAML to archive workflows.")
        # importing requests locally because it interferes with instantiating conda environments
        import requests

        env_archive = self.get_archive(dag)
        if os.path.exists(env_archive):
            return env_archive

        try:
            # Download
            logger.info("Downloading packages for conda environment {}...".format(self.file))
            os.makedirs(env_archive, exist_ok=True)
            try:
                out = subprocess.check_output(["conda", "list", "--explicit",
                    "--prefix", self.get_path(dag)],
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

    def get_path(self, dag):
        """Return environment path from hash.
        First tries full hash, if it does not exist, (8-prefix) is used as
        default."""
        for h in [self.hash, self.hash[:8]]:
            path = os.path.join(dag.workflow.persistence.conda_env_path, h)
            if os.path.exists(path):
                return path
        return path

    def create(self, dag):
        """ Create the conda enviroment."""
        # Read env file and create hash.
        env_file = self.file
        tmp_file = None

        url_scheme, *_ = urlparse(env_file)
        if url_scheme and not url_scheme == 'file':
            with tempfile.NamedTemporaryFile(delete=False) as tmp:
                tmp.write(self.content)
                env_file = tmp.name
                tmp_file = tmp.name

        env_hash = self.hash
        env_path = self.get_path(dag)
        # Create environment if not already present.
        if not os.path.exists(env_path):
            logger.info("Creating conda environment {}...".format(env_file))
            # Check if env archive exists. Use that if present.
            env_archive = self.get_archive(dag)
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
                logger.info("Environment for {} created.".format(env_file))
            except subprocess.CalledProcessError as e:
                # remove potential partially installed environment
                shutil.rmtree(env_path, ignore_errors=True)
                raise CreateCondaEnvironmentException(
                    "Could not create conda environment from {}:\n".format(env_file) +
                    e.output.decode())

        if tmp_file:
            # temporary file was created
            os.remove(tmp_file)

        self.path = env_path
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
