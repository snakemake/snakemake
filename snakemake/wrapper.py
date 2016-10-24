__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os

from snakemake.script import script


def is_script(path):
    return path.endswith("wrapper.py") or path.endswith("wrapper.R")


def get_path(path):
    if not (path.startswith("http") or path.startswith("file:")):
        path = "https://bitbucket.org/snakemake/snakemake-wrappers/raw/" + path
    return path


def get_script(path):
    path = get_path(path)
    if not is_script(path):
        path += "/wrapper.py"
    return path


def get_conda_env(path):
    path = get_path(path)
    if is_script(path):
        # TODO rather use dirname equivalent for urls
        path = os.path.dirname(path)
    return path + "environment.yaml"


def wrapper(path, input, output, params, wildcards, threads, resources, log, config, conda_env):
    """
    Load a wrapper from https://bitbucket.org/snakemake/snakemake-wrappers under
    the given path + wrapper.py and execute it.
    """
    path = get_script(path)
    script(path, "", input, output, params, wildcards, threads, resources, log, config, conda_env)
