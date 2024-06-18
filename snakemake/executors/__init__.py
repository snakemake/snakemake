__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import contextlib

from snakemake.logging import logger


@contextlib.contextmanager
def change_working_directory(directory=None):
    """Change working directory in execution context if provided."""
    if directory:
        try:
            saved_directory = os.getcwd()
            logger.info(f"Changing to shadow directory: {directory}")
            os.chdir(directory)
            yield
        finally:
            os.chdir(saved_directory)
    else:
        yield
