__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


from snakemake.exceptions import WorkflowError
from snakemake.script import script
from snakemake.sourcecache import SourceCache, infer_source_file


PREFIX = "https://github.com/snakemake/snakemake-wrappers/raw/"

EXTENSIONS = [".py", ".R", ".Rmd", ".jl"]


def is_script(source_file):
    filename = source_file.get_filename()
    return (
        filename.endswith("wrapper.py")
        or filename.endswith("wrapper.R")
        or filename.endswith("wrapper.jl")
    )


def get_path(path, prefix=None):
    if not is_url(path):
        if prefix is None:
            prefix = PREFIX
        elif prefix.startswith("git+file"):
            parts = path.split("/")
            path = "/" + "/".join(parts[1:]) + "@" + parts[0]
        path = prefix + path
    return infer_source_file(path)


def is_url(path):
    return (
        path.startswith("http")
        or path.startswith("file:")
        or path.startswith("git+file")
    )


def find_extension(source_file, sourcecache: SourceCache):
    for ext in EXTENSIONS:
        if source_file.get_filename().endswith(f"wrapper{ext}"):
            return source_file

    for ext in EXTENSIONS:
        script = source_file.join(f"wrapper{ext}")

        if sourcecache.exists(script):
            return script


def get_script(path, sourcecache: SourceCache, prefix=None):
    path = get_path(path, prefix=prefix)
    return find_extension(path, sourcecache)


def get_conda_env(path, prefix=None):
    path = get_path(path, prefix=prefix)
    if is_script(path):
        # URLs and posixpaths share the same separator. Hence use posixpath here.
        path = path.get_basedir()
    return path.join("environment.yaml")


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
    sourcecache_path,
    runtime_sourcecache_path,
):
    """
    Load a wrapper from https://github.com/snakemake/snakemake-wrappers under
    the given path + wrapper.(py|R|Rmd) and execute it.
    """
    assert path is not None
    script_source = get_script(
        path,
        SourceCache(sourcecache_path, runtime_cache_path=runtime_sourcecache_path),
        prefix=prefix,
    )
    if script_source is None:
        raise WorkflowError(
            f"Unable to locate wrapper script for wrapper {path}. "
            "This can be a network issue or a mistake in the wrapper URL."
        )
    script(
        script_source.get_path_or_uri(),
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
        sourcecache_path,
        runtime_sourcecache_path,
    )
