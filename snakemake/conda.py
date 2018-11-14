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

from snakemake.exceptions import CreateCondaEnvironmentException, WorkflowError
from snakemake.logging import logger
from snakemake.common import strip_prefix
from snakemake import utils
from snakemake import singularity
from snakemake.io import git_content


def content(env_file):
    if env_file.startswith("git+file:"):
        return git_content(env_file).encode('utf-8')
    elif urlparse(env_file).scheme:
        try:
            return urlopen(env_file).read()
        except URLError as e:
            raise WorkflowError("Failed to open environment file {}:".format(env_file), e)
    else:
        if not os.path.exists(env_file):
            raise WorkflowError("Conda env file does not "
                                "exist: {}".format(env_file))
        with open(env_file, 'rb') as f:
            return f.read()


class Env:

    """Conda environment from a given specification file."""

    def __init__(self, env_file, dag, singularity_img=None):
        self.file = env_file

        self._env_dir = dag.workflow.persistence.conda_env_path
        self._env_archive_dir = dag.workflow.persistence.conda_env_archive_path

        self._hash = None
        self._content_hash = None
        self._content = None
        self._path = None
        self._archive_file = None
        self._singularity_img = singularity_img

    @property
    def singularity_img_url(self):
        return self._singularity_img.url if self._singularity_img else None

    @property
    def content(self):
        if self._content is None:
            self._content = content(self.file)
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
            if self._singularity_img:
                md5hash.update(self._singularity_img.url.encode())
            md5hash.update(self.content)
            self._hash = md5hash.hexdigest()
        return self._hash

    @property
    def content_hash(self):
        if self._content_hash is None:
            md5hash = hashlib.md5()
            md5hash.update(self.content)
            self._content_hash = md5hash.hexdigest()
        return self._content_hash

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
            self._archive_file = os.path.join(self._env_archive_dir, self.content_hash)
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
            with open(os.path.join(env_archive,
                                   "packages.txt"), "w") as pkg_list:
                for l in out.decode().split("\n"):
                    if l and not l.startswith("#") and not l.startswith("@"):
                        pkg_url = l
                        logger.info(pkg_url)
                        parsed = urlparse(pkg_url)
                        pkg_name = os.path.basename(parsed.path)
                        # write package name to list
                        print(pkg_name, file=pkg_list)
                        # download package
                        pkg_path = os.path.join(env_archive, pkg_name)
                        with open(pkg_path, "wb") as copy:
                            r = requests.get(pkg_url)
                            r.raise_for_status()
                            copy.write(r.content)
                        try:
                            tarfile.open(pkg_path)
                        except:
                            raise WorkflowError("Package is invalid tar archive: {}".format(pkg_url))
        except (requests.exceptions.ChunkedEncodingError, requests.exceptions.HTTPError) as e:
            shutil.rmtree(env_archive)
            raise WorkflowError("Error downloading conda package {}.".format(pkg_url))
        except (Exception, BaseException) as e:
            shutil.rmtree(env_archive)
            raise e
        return env_archive

    def create(self, dryrun=False):
        """ Create the conda enviroment."""
        if self._singularity_img:
            check_conda(self._singularity_img)

        # Read env file and create hash.
        env_file = self.file
        tmp_file = None

        url_scheme, *_ = urlparse(env_file)
        if (url_scheme and not url_scheme == 'file') or \
            (not url_scheme and env_file.startswith("git+file:/")):
            with tempfile.NamedTemporaryFile(delete=False, suffix=".yaml") as tmp:
                tmp.write(self.content)
                env_file = tmp.name
                tmp_file = tmp.name

        env_hash = self.hash
        env_path = self.path

        # Check for broken environment
        if os.path.exists(os.path.join(env_path,"env_setup_start")) and not os.path.exists(os.path.join(env_path,"env_setup_done")):
            if dryrun:
                logger.info("Incomplete Conda environment {} will be recreated.".format(utils.simplify_path(self.file)))
            else:
                logger.info("Removing incomplete Conda environment {}...".format(utils.simplify_path(self.file)))
                shutil.rmtree(env_path, ignore_errors=True)

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
                # Touch "start" flag file
                os.makedirs(env_path, exist_ok=True)
                with open(os.path.join(env_path,"env_setup_start"), "a") as f:
                    pass

                if os.path.exists(env_archive):
                    logger.info("Using archived local conda packages.")
                    pkg_list = os.path.join(env_archive, "packages.txt")
                    if os.path.exists(pkg_list):
                        # read pacakges in correct order
                        # this is for newer env archives where the package list
                        # was stored
                        packages = [os.path.join(env_archive, pkg.rstrip())
                                    for pkg in open(pkg_list)]
                    else:
                        # guess order
                        packages = glob(os.path.join(env_archive, "*.tar.bz2"))

                    # install packages manually from env archive
                    cmd = " ".join(
                        ["conda", "create", "--copy", "--prefix", env_path] +
                        packages)
                    if self._singularity_img:
                        cmd = singularity.shellcmd(self._singularity_img.path, cmd)
                    out = subprocess.check_output(cmd, shell=True,
                                                  stderr=subprocess.STDOUT)

                else:
                    # Copy env file to env_path (because they can be on
                    # different volumes and singularity should only mount one).
                    # In addition, this allows to immediately see what an
                    # environment in .snakemake/conda contains.
                    target_env_file = env_path + ".yaml"
                    shutil.copy(env_file, target_env_file)

                    logger.info("Downloading remote packages.")
                    cmd = " ".join(["conda", "env", "create",
                                                "--file", target_env_file,
                                                "--prefix", env_path])
                    if self._singularity_img:
                        cmd = singularity.shellcmd(self._singularity_img.path, cmd)
                    out = subprocess.check_output(cmd, shell=True,
                                                  stderr=subprocess.STDOUT)
                # Touch "done" flag file
                with open(os.path.join(env_path,"env_setup_done"), "a") as f:
                    pass

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

    def __hash__(self):
        # this hash is only for object comparison, not for env paths
        return hash(self.file)

    def __eq__(self, other):
        if isinstance(other, Env):
            return self.file == other.file
        return False


def shellcmd(env_path):
    return "source activate {};".format(env_path)


def check_conda(singularity_img=None):
    def get_cmd(cmd):
        if singularity_img:
            return singularity.shellcmd(singularity_img.path, cmd)
        return cmd

    try:
        # Use type here since conda now is a function.
        # type allows to check for both functions and regular commands.
        subprocess.check_output(get_cmd("type conda"),
                                shell=True,
                                stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        if singularity_img:
            raise CreateCondaEnvironmentException("The 'conda' command is not "
                                                  "available inside "
                                                  "your singularity container "
                                                  "image.")
        else:
            raise CreateCondaEnvironmentException("The 'conda' command is not "
                                                  "available.")
    try:
        version = subprocess.check_output(get_cmd("conda --version"),
                                          shell=True,
                                          stderr=subprocess.STDOUT).decode() \
                                                                   .split()[1]
        if StrictVersion(version) < StrictVersion("4.2"):
            raise CreateCondaEnvironmentException(
                "Conda must be version 4.2 or later."
            )
    except subprocess.CalledProcessError as e:
        raise CreateCondaEnvironmentException(
            "Unable to check conda version:\n" + e.output.decode()
        )
