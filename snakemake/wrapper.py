__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
import posixpath

from urllib.error import URLError
from urllib.request import urlopen

from snakemake.exceptions import WorkflowError
from snakemake.script import script

# where to look for type of script?

def is_script(path):
    return path.endswith("wrapper.py") or path.endswith("wrapper.R")


def get_path(path, prefix=None):
    if not (path.startswith("http") or path.startswith("file:")):
        if prefix is None:
            prefix = "https://bitbucket.org/snakemake/snakemake-wrappers/raw/"
        path = prefix + path
    return path

def find_local_extension(path, possible_extensions_ordered):
    for ext in possible_extensions_ordered:
        wrapper_with_ext = "/wrapper" + ext
        if os.path.exists(path.replace("file://", "") + wrapper_with_ext):
            path += wrapper_with_ext
            break
    else:
        exts = "/".join(possible_extensions_ordered)
        raise WorkflowError("No local wrapper {path}/wrapper({exts}) found".format(path=path, exts=exts))

    return path

def find_remote_extension(path, possible_extensions_ordered):
    for ext in possible_extensions_ordered:
        wrapper_with_ext = "/wrapper" + ext
        try:
            urlopen(path + wrapper_with_ext)
            path += wrapper_with_ext
            break
        except URLError:
            continue
    else:
        exts = "/".join(possible_extensions_ordered)
        raise WorkflowError("No remote wrapper {path}/wrapper({exts}) found".format(path=path, exts=exts))

    return path

def get_script(path, prefix=None):
    path = get_path(path, prefix=prefix)

    is_local = path.startswith("file:")
    possible_extensions_ordered = [".py", ".R", ".Rmd"]

    # chesk for path locally
    if not is_script(path) and is_local:
        path = find_local_extension(path, possible_extensions_ordered)
    elif not is_script(path) and not is_local:
        path = find_remote_extension(path, possible_extensions_ordered)

    return path


def get_conda_env(path, prefix=None):
    path = get_path(path, prefix=prefix)
    if is_script(path):
        # URLs and posixpaths share the same separator. Hence use posixpath here.
        path = posixpath.dirname(path)
    return path + "/environment.yaml"


def wrapper(path,
            input,
            output,
            params,
            wildcards,
            threads,
            resources,
            log,
            config,
            rulename,
            conda_env,
            singularity_img,
            singularity_args,
            bench_record,
            prefix):
    """
    Load a wrapper from https://bitbucket.org/snakemake/snakemake-wrappers under
    the given path + wrapper.py and execute it.
    """
    path = get_script(path, prefix=prefix)
    script(path, "", input, output, params, wildcards, threads, resources,
           log, config, rulename, conda_env, singularity_img,
           singularity_args, bench_record)
