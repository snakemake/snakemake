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
from snakemake.sourcecache import SourceCache, infer_source_file


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


def find_extension(
    path, sourcecache: SourceCache, extensions=[".py", ".R", ".Rmd", ".jl"]
):
    for ext in extensions:
        if path.endswith("wrapper{}".format(ext)):
            return path

    path = infer_source_file(path)
    for ext in extensions:
        script = path.join("wrapper{}".format(ext))

        if sourcecache.exists(script):
            return script


def get_script(path, sourcecache: SourceCache, prefix=None):
    path = get_path(path, prefix=prefix)
    return find_extension(path, sourcecache)


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
    conda_base_path,
    container_img,
    singularity_args,
    env_modules,
    bench_record,
    prefix,
    jobid,
    bench_iteration,
    cleanup_scripts,
    shadow_dir,
    runtime_sourcecache_path,
):
    """
    Load a wrapper from https://github.com/snakemake/snakemake-wrappers under
    the given path + wrapper.(py|R|Rmd) and execute it.
    """
    path = get_script(
        path, SourceCache(runtime_cache_path=runtime_sourcecache_path), prefix=prefix
    )
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
        conda_base_path,
        container_img,
        singularity_args,
        env_modules,
        bench_record,
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
        runtime_sourcecache_path,
    )
