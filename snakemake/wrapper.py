__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


import os
import posixpath

from urllib.error import URLError
from urllib.request import urlopen

from snakemake.exceptions import WorkflowError
from snakemake.script import script


PREFIX = "https://github.com/snakemake/snakemake-wrappers/raw/"


def is_script(path):
    return (
        path.endswith("wrapper.py")
        or path.endswith("wrapper.R")
        or path.endswith("wrapper.jl")
    )


def get_path(path, prefix=None):
    if not is_url(path):
        if prefix is None:
            prefix = PREFIX
        elif prefix.startswith("git+file"):
            parts = path.split("/")
            path = "/" + "/".join(parts[1:]) + "@" + parts[0]
        path = prefix + path
    return path


def is_url(path):
    return (
        path.startswith("http")
        or path.startswith("file:")
        or path.startswith("git+file")
    )


def is_local(path):
    return path.startswith("file:")


def is_git_path(path):
    return path.startswith("git+file:")


def find_extension(path, extensions=[".py", ".R", ".Rmd", ".jl"]):
    for ext in extensions:
        if path.endswith("wrapper{}".format(ext)):
            return path
    for ext in extensions:
        script = "/wrapper{}".format(ext)
        if is_local(path):
            if path.startswith("file://"):
                p = path[7:]
            elif path.startswith("file:"):
                p = path[5:]
            if os.path.exists(p + script):
                return path + script
        else:
            try:
                urlopen(path + script)
                return path + script
            except URLError:
                continue
    if is_git_path(path):
        path, version = path.split("@")
        return os.path.join(path, "wrapper.py") + "@" + version
    else:
        return path + "/wrapper.py"  # default case


def get_script(path, prefix=None):
    path = get_path(path, prefix=prefix)
    return find_extension(path)


def get_conda_env(path, prefix=None):
    path = get_path(path, prefix=prefix)
    if is_script(path):
        # URLs and posixpaths share the same separator. Hence use posixpath here.
        path = posixpath.dirname(path)
    if is_git_path(path):
        path, version = path.split("@")
        return os.path.join(path, "environment.yaml") + "@" + version
    return path + "/environment.yaml"


def wrapper(
    path,
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
    container_img,
    singularity_args,
    env_modules,
    bench_record,
    prefix,
    jobid,
    bench_iteration,
    cleanup_scripts,
    shadow_dir,
):
    """
    Load a wrapper from https://github.com/snakemake/snakemake-wrappers under
    the given path + wrapper.(py|R|Rmd) and execute it.
    """
    path = get_script(path, prefix=prefix)
    script(
        path,
        "",
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
        container_img,
        singularity_args,
        env_modules,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
    )
