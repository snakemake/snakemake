from typing import TYPE_CHECKING

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import re
from typing import Optional, Dict
from snakemake.exceptions import WorkflowError
from snakemake.script import script
from snakemake.sourcecache import GithubFile, SourceCache, infer_source_file
from snakemake.executors.local import RunArgs
from snakemake.sourcecache import SourceFile
from snakemake_interface_software_deployment_plugins import EnvSpecSourceFile

if TYPE_CHECKING:
    from snakemake_software_deployment_plugin_conda import EnvSpec as CondaEnvSpec

EXTENSIONS = [".py", ".R", ".Rmd", ".jl"]

DEFAULT_WRAPPER_PREFIX = "https://github.com/snakemake/snakemake-wrappers/raw/"


ver_regex = re.compile(r"v?(?P<ver>[0-9]+\.[0-9]+\.[0-9]+)")


def is_script(source_file):
    filename = source_file.get_filename()
    return (
        filename.endswith("wrapper.py")
        or filename.endswith("wrapper.R")
        or filename.endswith("wrapper.jl")
    )


def get_path(path: str, prefix: Optional[str] = None):
    if not is_url(path):
        if prefix is not None and prefix.startswith("git+file"):
            parts = path.split("/")
            path = "/" + "/".join(parts[1:]) + "@" + parts[0]
        if prefix is None:
            ref, path_tail = path.split("/", 1)
            if ver_regex.match(ref):
                # Best case, we have a github file with tag, so that we can
                # persistently cache it in the sourcecache.
                return GithubFile(
                    repo="snakemake/snakemake-wrappers", tag=ref, path=path_tail
                )
            else:
                # Otherwise, use a plain url and store it in runtime cache only.
                path = DEFAULT_WRAPPER_PREFIX + path
        else:
            path = prefix + path
    return infer_source_file(path)


def is_url(path):
    return (
        path.startswith("http")
        or path.startswith("file:")
        or path.startswith("git+file")
    )


def find_extension(
    source_file, sourcecache: SourceCache
) -> SourceFile | Dict[str, str]:
    for ext in EXTENSIONS:
        if source_file.get_filename().endswith(f"wrapper{ext}"):
            return source_file

    errors = {}
    for ext in EXTENSIONS:
        script_name = f"wrapper{ext}"
        script = source_file.join(script_name)

        try:
            sourcecache.try_access(script)
            return script
        except Exception as e:
            if isinstance(e, WorkflowError):
                msg = str(e)
            else:
                msg = e.__class__.__name__
                if str(e):
                    msg += f": {str(e)}"
            errors[script_name] = msg
    return errors


def get_script(
    path, sourcecache: SourceCache, prefix=None
) -> SourceFile | Dict[str, str]:
    path = get_path(path, prefix=prefix)
    return find_extension(path, sourcecache)


def get_conda_env(path, sourcecache: SourceCache, prefix=None) -> "CondaEnvSpec":
    from snakemake_software_deployment_plugin_conda import EnvSpec as CondaEnvSpec

    path = get_path(path, prefix=prefix)
    if is_script(path):
        # URLs and posixpaths share the same separator. Hence use posixpath here.
        path = path.get_basedir()
    path = path.join("environment.yaml")
    cached_path = sourcecache.cache(path)
    spec = CondaEnvSpec(
        envfile=EnvSpecSourceFile(
            path_or_uri=path.get_path_or_uri(secret_free=True), cached=cached_path
        )
    )
    spec.technical_init()
    return spec


def wrapper(
    path,
    run_args: RunArgs,
    config: Dict,
) -> None:
    """
    Load a wrapper from https://github.com/snakemake/snakemake-wrappers under
    the given path + wrapper.(py|R|Rmd) and execute it.
    """
    assert path is not None
    script_source = get_script(
        path,
        SourceCache(
            run_args.cache_path, runtime_cache_path=run_args.runtime_cache_path
        ),
        prefix=run_args.wrapper_prefix,
    )
    if not isinstance(script_source, SourceFile):
        raise WorkflowError(
            f"Unable to locate wrapper script at {run_args.wrapper_prefix}{path}. "
            "This can be a network issue or a mistake in the wrapper URL. Encountered errors:\n"
            + "\n".join(f"{script}: {error}" for script, error in script_source.items())
        )
    script(
        script_source,
        run_args,
        config,
    )
