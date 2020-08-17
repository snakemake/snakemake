__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016-2020, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
import posixpath
import re
import hashlib
import requests
import yaml

from urllib.error import URLError
from urllib.request import urlopen

from snakemake.deployment.docker import DockerContainer
from snakemake.exceptions import WorkflowError
from snakemake.script import script

# Derive custom container prefix from the environment, default to quay.io
CONTAINER_PREFIX = os.environ.get(
    "SNAKEMAKE_WRAPPER_CONTAINER_PREFIX", "quay.io/snakemake-wrappers"
)
WRAPPER_REPOSITORY = os.environ.get(
    "SNAKEMAKE_WRAPPER_REPOSITORY",
    "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/",
)

# File prefix is needed to eliminate code smell
file_prefix = "file:"


class WrapperSpec:
    """A Snakemake Wrapper spec coincides with a snakemake wrapper. Specifically,
       we can use it to interact with scripts (local or on GitHub) and can
       generate a unique hash and container name based on the environment.yaml.
       We can then easily determine if the container exists.
    """

    def __init__(self, path, prefix=None):
        self.prefix = prefix

        # Populated with properties
        self._meta = None
        self._environment = None
        self._environment_file = None
        self._packages = {}
        self._path = path

    @property
    def conda_env(self):
        """Get the conda environment
        """
        if self._environment_file is None:
            path = self.path
            if is_git_path(path):
                path, version = path.split("@")
                self._environment_file = (
                    os.path.join(path, "environment.yaml") + "@" + version
                )
            else:
                self._environment_file = path + "/environment.yaml"
        return self._environment_file

    @property
    def packages(self):
        """Read the full environment.yaml file
        """
        environment_file = self.conda_env
        if not environment_file:
            return

        # Read the file, whether local or remote, and load yaml
        self._environment = self.read_local_remote_file(environment_file)
        self._packages = yaml.load(self._environment, Loader=yaml.SafeLoader)
        return self._packages

    @property
    def path(self):
        prefix = self.prefix or WRAPPER_REPOSITORY
        path = get_path(self._path, prefix=prefix)
        if is_script(path):
            path = posixpath.dirname(path)
        return path

    def parse_path(self, path):
        """Given a path, determine if it has a script associated. We set 
           self._path to be the directory (web or local) name, and self._script
           to be the wrapper.
        """
        # A script can be online (https) or local
        if is_script(path):
            self._script = path.split("/")[-1]
            self._path = "/".join(path.split("/")[:-1])
        else:
            self._path = path

    @property
    def meta(self):
        """Given a path to a rule, load the meta.yaml / meta.yml and determine if
           a container is provided. If so, return the container URI. Otherwise
           return None. 
        """
        if self._meta is not None:
            return self._meta

        if is_git_path(self.path):
            path, version = self.path.split("@")
            path = os.path.join(path, "meta.yaml") + "@" + version
        else:
            path = os.path.join(self.path, "meta.yaml")

        # Read the file, whether local or remote
        meta = self.read_local_remote_file(path)

        # If meta doesn't exist (local or remote) we return None
        if meta:
            self._meta = yaml.load(meta, Loader=yaml.SafeLoader)
        return self._meta

    @property
    def sets(self):
        """Return sorted package sets (name and version)
        """
        # Generate hashes for the package names (container name) and tag (versions)
        # {'bwa': '0.7.17', 'samtools': '1.9', 'picard': '2.20.1'}
        sets = {
            d[0].strip(): d[1].strip()
            for d in [
                re.split("==|=|<=|>=", x) for x in self.packages.get("dependencies", [])
            ]
        }
        return sorted(sets.items())

    def read_local_remote_file(self, path):
        """Given a filename (local or remote) read it with a web request
           or directly from the filesystem.
        """
        if path.startswith(file_prefix):
            path = re.sub("(file:|file:/|file://)", "", path).strip("/")

        # If the file exists, read from filesystem
        content = None
        if os.path.exists(path):
            with open(path, "r") as fd:
                content = fd.read()

        # Otherwise, if it's a url, retrieve and read it
        elif path.startswith("http"):
            try:
                content = urlopen(path).read()
            except URLError:
                return
        return content

    def _get_md5(self, listing):
        """Private function to generate an md5sum for a listing of variables
           You should use the get_container_shamd5, get_versions_md5, and
           get_names_md5 for more finished hashes. 
        """
        hasher = hashlib.md5()
        hasher.update("".join(listing).encode("utf-8"))
        return hasher.hexdigest()

    def get_container_name(self):
        """Given a base container name, and hashes for packages and names,
           return the full container unique resource identifier.
        """
        container_name = None
        if self.meta is not None:
            spec = self.meta.get("container", {}).get("spec")

            # Case 1: meta specifies that the container builds to snakemake-wrappers (hash)
            if spec == "auto":
                container_name = "%s/%s:%s" % (
                    CONTAINER_PREFIX,
                    self.get_names_md5(),
                    self.get_versions_md5(),
                )
            # Case 2: developer wants to link to pre-existing container
            elif spec is not None and re.search("docker:", spec):
                container_name = spec.replace("docker://", "", 1)

            # Case 3 there is a Dockerfile - I'm not sure how we would get container name from it (and need to build it?)
            # Name is None if container Doesn't exist
            return DockerContainer(container_name).uri

    def get_versions_md5(self):
        """Given packages for a wrapper, return the md5 sum of sorted versions
           (sorted by package names)
        """
        # Index at 1 is the version
        return self._get_md5([x[1] for x in self.sets])

    def get_names_md5(self):
        # Index at 0 is the name
        return self._get_md5([x[0] for x in self.sets])


def get_script(path, prefix=None):
    """Get script uses the original (unparsed) path.
    """
    path = get_path(path, prefix=prefix)
    return find_extension(path)


def is_script(path):
    return (
        path.endswith("wrapper.py")
        or path.endswith("wrapper.R")
        or path.endswith("wrapper.jl")
    )


def get_path(path, prefix=None):
    if not (
        path.startswith("http")
        or path.startswith(file_prefix)
        or path.startswith("git+file")
    ):
        if prefix is None:
            prefix = "https://github.com/snakemake/snakemake-wrappers/raw/"
        elif prefix.startswith("git+file"):
            parts = path.split("/")
            path = "/" + "/".join(parts[1:]) + "@" + parts[0]
        path = prefix + path
    return path


def is_local(path):
    return path.startswith(file_prefix)


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
            elif path.startswith(file_prefix):
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
