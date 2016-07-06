import os
import shutil
import subprocess
import json
import tempfile

from snakemake.exceptions import CreateEnvironmentException
from snakemake.logging import logger


class Environments:

    def __init__(self):
        self.path = os.path.abspath(".conda")
        self.environments = dict()

    def register(self, env_file):
        if env_file not in self.environments:
            self.environments[env_file] = tempfile.mkdtemp(dir=self.path)

    def cleanup(self):
        for _, envdir in self.environments.items():
            shutil.rmtree(env, ignore_errors=True)

    def create(self, env_file):
        """Create conda environment if specified."""
        self.register(env_file)
        env = self[env_file]
        temp_env_file = None
        if not os.path.exists(env):
            logger.info("Creating conda environment for {}...".format(env_file))
            os.makedirs(os.path.dirname(env), exist_ok=True)
            try:
                remote = urlopen(env_file)
                with tempfile.NamedTemporaryFile(delete=False) as temp_env_file:
                    temp_env_file.write(remote.read())
                    env_file = temp_env_file.name
            except ValueError:
                # no download necessary
                pass
            try:
                out = subprocess.check_output(["conda", "env", "create",
                                               "--file", env_file,
                                               "--prefix", env],
                                               stderr=subprocess.STDOUT)
                logger.info("Environment for {} created.".format(env_file))
            except subprocess.CalledProcessError as e:
                raise CreateEnvironmentException(
                    "Could not create conda environment from {}:\n".format(env_file) +
                    e.output.decode())
            if temp_env_file is not None:
                os.remove(env_file)

    def __getitem__(self, env_file):
        return self.environments[env_file]
