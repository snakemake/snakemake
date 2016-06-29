import os
import shutil
import subprocess
import json

from snakemake.exceptions import CreateEnvironmentException
from snakemake.logging import logger


class Environments:

    def __init__(self):
        self.path = os.path.abspath(".conda")
        self.environments = dict()

    def register(self, env_file):
        env = os.path.abspath(os.path.join(self.path, env_file))
        self.environments[env_file] = env
        if os.path.exists(env):
            shutil.rmtree(env)

    def create(self, env_file):
        """Create conda environment if specified."""
        if env_file not in self.environments:
            self.register(env_file)
        env = self[env_file]
        if not os.path.exists(env):
            logger.info("Creating conda environment for {}...".format(env_file))
            os.makedirs(os.path.dirname(env), exist_ok=True)
            try:
                out = subprocess.check_output(["conda", "env", "create",
                                               "--file", env_file,
                                               "--prefix", env],
                                               stderr=subprocess.STDOUT)
                logger.info("Environment for {} created.".format(env_file))
            except subprocess.CalledProcessError as e:
                raise CreateEnvironmentException(
                    "Could not create conda environment from {}:\n".format(env_file) +
                    json.loads(e.output.decode())["error"])

    def __getitem__(self, env_file):
        return self.environments[env_file]
