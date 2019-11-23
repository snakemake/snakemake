__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
from urllib.request import urlopen, pathname2url

from snakemake.utils import format
from snakemake.io import git_content, split_git_path
from snakemake.script_handlers import get_executor_class


def get_source(path, basedir="."):
    source = None
    if not path.startswith("http") and not path.startswith("git+file"):
        if path.startswith("file://"):
            path = path[7:]
        elif path.startswith("file:"):
            path = path[5:]
        if not os.path.isabs(path):
            path = os.path.abspath(os.path.join(basedir, path))
        path = "file://" + path
    path = format(path, stepout=1)
    if path.startswith("file://"):
        sourceurl = "file:" + pathname2url(path[7:])
    elif path.startswith("git+file"):
        source = git_content(path)
        (root_path, file_path, version) = split_git_path(path)
        path = path.rstrip("@" + version)
    else:
        sourceurl = path

    language = None
    if path.endswith(".py"):
        language = "python"
    elif path.endswith(".ipynb"):
        language = "jupyter"
    elif path.endswith(".R"):
        language = "r"
    elif path.endswith(".Rmd"):
        language = "rmarkdown"
    elif path.endswith(".jl"):
        language = "julia"

    if source is None:
        with urlopen(sourceurl) as source:
            return path, source.read(), language
    else:
        return path, source, language


def script(
    path,
    basedir,
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
    jobid,
    bench_iteration,
    cleanup_scripts,
    shadow_dir,
):
    """
    Load a script from the given basedir + path and execute it.
    """
    path, source, language = get_source(path, basedir)

    ExecClass = get_executor_class(language)
    executor = ExecClass(
        path,
        source,
        basedir,
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
        jobid,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
    )

    executor.evaluate()
