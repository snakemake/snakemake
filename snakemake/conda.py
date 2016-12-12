import os
import subprocess
import tempfile
from urllib.request import urlopen
import hashlib
import shutil
from distutils.version import StrictVersion

from snakemake.exceptions import CreateCondaEnvironmentException
from snakemake.logging import logger


def create_env(job):
    """ Create conda enviroment for the given job. """
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

    md5hash = hashlib.md5()
    env_file = job.conda_env_file
    if os.path.exists(env_file):
        with open(env_file, 'rb') as f:
            md5hash.update(f.read())
    else:
        if env_file.startswith('file://'):
            with open(env_file.replace('file://', ''), 'rb') as handle:
                    content = handle.read()
        else:
            content = urlopen(env_file).read()
        md5hash.update(content)
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp.write(content)
            env_file = tmp.name

    env_path = os.path.join(job.rule.workflow.persistence.conda_env_path, md5hash.hexdigest()[:8])
    if not os.path.exists(env_path):
        logger.info("Creating conda environment for {}...".format(job.conda_env_file))
        try:
            out = subprocess.check_output(["conda", "env", "create",
                                           "--file", env_file,
                                           "--prefix", env_path],
                                           stderr=subprocess.STDOUT)
            logger.debug(out)
            logger.info("Environment for {} created.".format(job.conda_env_file))
        except subprocess.CalledProcessError as e:
            # remove potential partially installed environment
            shutil.rmtree(env_path, ignore_errors=True)
            raise CreateCondaEnvironmentException(
                "Could not create conda environment from {}:\n".format(job.conda_env_file) +
                e.output.decode())

    if env_file != job.conda_env_file:
        # temporary file was created
        os.remove(env_file)

    return env_path
