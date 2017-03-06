import os
import subprocess
import tempfile
from urllib.request import urlopen, urlretrieve
from urllib.parse import urlparse
import hashlib
import shutil
from distutils.version import StrictVersion
import json
from glob import glob

from snakemake.exceptions import CreateCondaEnvironmentException, WorkflowError
from snakemake.logging import logger


def get_env_archive(job, env_hash):
    """Get path to archived environment derived from given environment file."""
    return os.path.join(job.rule.workflow.persistence.conda_env_archive_path, env_hash)


def archive_env(job):
    """Create self-contained archive of environment."""
    try:
        import yaml
    except ImportError:
        raise WorkflowError("Error importing PyYAML. "
            "Please install PyYAML to archive workflows.")
    # importing requests locally because it interferes with instantiating conda environments
    import requests

    env_archive = get_env_archive(job, get_env_hash(job.conda_env_file))
    if os.path.exists(env_archive):
        return env_archive

    try:
        # Download
        logger.info("Downloading packages for conda environment {}...".format(job.conda_env_file))
        os.makedirs(env_archive, exist_ok=True)
        try:
            out = subprocess.check_output(["conda", "list", "--explicit",
                "--prefix", job.conda_env],
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


def is_remote_env_file(env_file):
    return urlparse(env_file).scheme


def get_env_hash(env_file):
    md5hash = hashlib.md5()
    if is_remote_env_file(env_file):
        md5hash.update(urlopen(env_file).read())
    else:
        with open(env_file, 'rb') as f:
            md5hash.update(f.read())
    return md5hash.hexdigest()


def get_env_path(job, env_hash):
    """Return environment path from hash.
    First tries full hash, if it does not exist, (8-prefix) is used as
    default."""
    for h in [env_hash, env_hash[:8]]:
        path = os.path.join(job.rule.workflow.persistence.conda_env_path, h)
        if os.path.exists(path):
            return path
    return path

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

def create_env(job):
    """ Create conda enviroment for the given job. """

    # Read env file and create hash.
    env_file = job.conda_env_file
    tmp_file = None
    is_remote = is_remote_env_file(env_file)
    if is_remote and is_remote != 'file':
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp.write(urlopen(env_file).read())
            env_file = tmp.name
            tmp_file = tmp.name
    env_hash = get_env_hash(env_file)
    env_path = get_env_path(job, env_hash)
    # Create environment if not already present.
    if not os.path.exists(env_path):
        logger.info("Creating conda environment {}...".format(job.conda_env_file))
        # Check if env archive exists. Use that if present.
        env_archive = get_env_archive(job, env_hash)
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
            logger.info("Environment for {} created.".format(job.conda_env_file))
        except subprocess.CalledProcessError as e:
            # remove potential partially installed environment
            shutil.rmtree(env_path, ignore_errors=True)
            raise CreateCondaEnvironmentException(
                "Could not create conda environment from {}:\n".format(job.conda_env_file) +
                e.output.decode())

    if tmp_file:
        # temporary file was created
        os.remove(tmp_file)

    return env_path
