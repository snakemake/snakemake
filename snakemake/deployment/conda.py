__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
from pathlib import Path
import re
from snakemake.sourcecache import (
    LocalGitFile,
    LocalSourceFile,
    SourceFile,
    infer_source_file,
)
import subprocess
import tempfile
from urllib.request import urlopen
from urllib.error import URLError
import hashlib
import shutil
from distutils.version import StrictVersion
import json
from glob import glob
import tarfile
import zipfile
import uuid
from enum import Enum
import threading
import shutil
from abc import ABC, abstractmethod


from snakemake.exceptions import CreateCondaEnvironmentException, WorkflowError
from snakemake.logging import logger
from snakemake.common import (
    is_local_file,
    lazy_property,
    parse_uri,
    strip_prefix,
    ON_WINDOWS,
)
from snakemake import utils
from snakemake.deployment import singularity, containerize
from snakemake.io import (
    IOFile,
    apply_wildcards,
    contains_wildcard,
    git_content,
    _IOFile,
)


class CondaCleanupMode(Enum):
    tarballs = "tarballs"
    cache = "cache"

    def __str__(self):
        return self.value


class Env:

    """Conda environment from a given specification file."""

    def __init__(
        self,
        workflow,
        env_file=None,
        env_name=None,
        env_dir=None,
        container_img=None,
        cleanup=None,
    ):
        self.file = env_file
        if env_file is not None:
            self.file = infer_source_file(env_file)
        self.name = env_name
        if env_name is not None:
            assert env_file is None, "bug: both env_file and env_name specified"

        self.frontend = workflow.conda_frontend
        self.workflow = workflow

        self._container_img = container_img
        self._env_dir = env_dir or (
            containerize.CONDA_ENV_PATH
            if self.is_containerized
            else workflow.persistence.conda_env_path
        )
        self._hash = None
        self._content_hash = None
        self._content = None
        self._content_deploy = None
        self._content_pin = None
        self._path = None
        self._archive_file = None
        self._cleanup = cleanup
        self._singularity_args = workflow.singularity_args

    @lazy_property
    def conda(self):
        return Conda(self._container_img)

    @lazy_property
    def pin_file(self):
        pin_file = Path(self.file.get_path_or_uri()).with_suffix(
            f".{self.conda.platform}.pin.txt"
        )

        if pin_file.exists():
            return infer_source_file(pin_file)
        else:
            return None

    @lazy_property
    def post_deploy_file(self):
        if self.file:
            deploy_file = Path(self.file.get_path_or_uri()).with_suffix(
                ".post-deploy.sh"
            )
            if deploy_file.exists():
                return infer_source_file(deploy_file)

    def _get_content(self):
        if self.is_named:
            from snakemake.shell import shell

            content = shell.check_output(
                "conda env export {}".format(self.address_argument),
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
            return content.encode()
        else:
            return self.workflow.sourcecache.open(self.file, "rb").read()

    def _get_content_deploy(self):
        self.check_is_file_based()
        if self.post_deploy_file:
            return self.workflow.sourcecache.open(self.post_deploy_file, "rb").read()
        return None

    def _get_content_pin(self):
        self.check_is_file_based()
        if self.pin_file:
            return self.workflow.sourcecache.open(self.pin_file, "rb").read()
        return None

    @property
    def _env_archive_dir(self):
        return self.workflow.persistence.conda_env_archive_path

    @property
    def container_img_url(self):
        return self._container_img.url if self._container_img else None

    @property
    def content(self):
        if self._content is None:
            self._content = self._get_content()
        return self._content

    @property
    def content_deploy(self):
        if self._content_deploy is None:
            self._content_deploy = self._get_content_deploy()
        return self._content_deploy

    @property
    def content_pin(self):
        if self._content_pin is None:
            self._content_pin = self._get_content_pin()
        return self._content_pin

    @property
    def hash(self):
        if self._hash is None:
            if self.is_containerized:
                self._hash = self.content_hash
            else:
                md5hash = hashlib.md5()
                # Include the absolute path of the target env dir into the hash.
                # By this, moving the working directory around automatically
                # invalidates all environments. This is necessary, because binaries
                # in conda environments can contain hardcoded absolute RPATHs.
                env_dir = os.path.realpath(self._env_dir)
                md5hash.update(env_dir.encode())
                if self._container_img:
                    md5hash.update(self._container_img.url.encode())
                content_deploy = self.content_deploy
                if content_deploy:
                    md5hash.update(content_deploy)
                md5hash.update(self.content)
                self._hash = md5hash.hexdigest()
        return self._hash

    @property
    def content_hash(self):
        if self._content_hash is None:
            md5hash = hashlib.md5()
            md5hash.update(self.content)
            content_deploy = self.content_deploy
            if content_deploy:
                md5hash.update(content_deploy)
            self._content_hash = md5hash.hexdigest()
        return self._content_hash

    @property
    def is_containerized(self):
        if not self._container_img:
            return False
        return self._container_img.is_containerized

    @property
    def is_named(self):
        return self.file is None

    def check_is_file_based(self):
        assert (
            self.file is not None
        ), "bug: trying to access conda env file based functionality for named environment"

    @property
    def address(self):
        """Path to directory of the conda environment.

        First tries full hash, if it does not exist, (8-prefix) is used
        as default.

        """
        if self.is_named:
            return self.name
        else:
            hash = self.hash
            env_dir = self._env_dir
            get_path = lambda h: os.path.join(env_dir, h)
            hash_candidates = [
                hash[:8],
                hash,
            ]  # [0] is the old fallback hash (shortened)
            exists = [os.path.exists(get_path(h)) for h in hash_candidates]
            if self.is_containerized or exists[1] or (not exists[0]):
                # containerizes, full hash exists or fallback hash does not exist: use full hash
                return get_path(hash_candidates[1])
            # use fallback hash
            return get_path(hash_candidates[0])

    @property
    def address_argument(self):
        if self.is_named:
            return "--name '{}'".format(self.address)
        else:
            return "--prefix '{}'".format(self.address)

    @property
    def archive_file(self):
        """Path to archive of the conda environment, which may or may not exist."""
        if self._archive_file is None:
            self._archive_file = os.path.join(self._env_archive_dir, self.content_hash)
        return self._archive_file

    def create_archive(self):
        """Create self-contained archive of environment."""
        from snakemake.shell import shell

        try:
            import yaml
        except ImportError:
            raise WorkflowError(
                "Error importing PyYAML. " "Please install PyYAML to archive workflows."
            )
        # importing requests locally because it interferes with instantiating conda environments
        import requests

        self.check_is_file_based()

        env_archive = self.archive_file
        if os.path.exists(env_archive):
            return env_archive

        try:
            # Download
            logger.info(
                "Downloading packages for conda environment {}...".format(
                    self.file.get_path_or_uri()
                )
            )
            os.makedirs(env_archive, exist_ok=True)
            try:
                out = shell.check_output(
                    "conda list --explicit {}".format(self.address_argument),
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                )
                logger.debug(out)
            except subprocess.CalledProcessError as e:
                raise WorkflowError("Error exporting conda packages:\n" + e.output)
            with open(os.path.join(env_archive, "packages.txt"), "w") as pkg_list:
                for l in out.split("\n"):
                    if l and not l.startswith("#") and not l.startswith("@"):
                        pkg_url = l
                        logger.info(pkg_url)
                        parsed = parse_uri(pkg_url)
                        pkg_name = os.path.basename(parsed.uri_path)
                        # write package name to list
                        print(pkg_name, file=pkg_list)
                        # download package
                        pkg_path = os.path.join(env_archive, pkg_name)
                        with open(pkg_path, "wb") as copy:
                            r = requests.get(pkg_url)
                            r.raise_for_status()
                            copy.write(r.content)
                        try:
                            if pkg_path.endswith(".conda"):
                                assert zipfile.ZipFile(pkg_path).testzip() is None
                            else:
                                tarfile.open(pkg_path)
                        except:
                            raise WorkflowError(
                                "Package is invalid tar/zip archive: {}".format(pkg_url)
                            )
        except (
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.HTTPError,
        ) as e:
            shutil.rmtree(env_archive)
            raise WorkflowError("Error downloading conda package {}.".format(pkg_url))
        except (Exception, BaseException) as e:
            shutil.rmtree(env_archive)
            raise e
        return env_archive

    def execute_deployment_script(self, env_file, deploy_file):
        """Execute post-deployment script if present"""
        from snakemake.shell import shell

        if ON_WINDOWS:
            raise WorkflowError(
                "Post deploy script {} provided for conda env {} but unsupported on windows.".format(
                    deploy_file, env_file
                )
            )
        logger.info(
            "Running post-deploy script {}...".format(
                os.path.relpath(path=deploy_file, start=os.getcwd())
            )
        )
        shell.check_output(
            self.conda.shellcmd(self.address, "sh {}".format(deploy_file)),
            stderr=subprocess.STDOUT,
        )

    def create(self, dryrun=False):
        """Create the conda enviroment."""
        from snakemake.shell import shell

        self.check_is_file_based()

        # Read env file and create hash.
        env_file = self.file
        deploy_file = None
        pin_file = None
        tmp_env_file = None
        tmp_deploy_file = None
        tmp_pin_file = None

        if not isinstance(env_file, LocalSourceFile) or isinstance(
            env_file, LocalGitFile
        ):
            with tempfile.NamedTemporaryFile(delete=False, suffix=".yaml") as tmp:
                # write to temp file such that conda can open it
                tmp.write(self.content)
                env_file = tmp.name
                tmp_env_file = tmp.name
            if self.post_deploy_file:
                with tempfile.NamedTemporaryFile(
                    delete=False, suffix=".post-deploy.sh"
                ) as tmp:
                    # write to temp file such that conda can open it
                    tmp.write(self.content_deploy)
                    deploy_file = tmp.name
                    tmp_deploy_file = tmp.name
            if self.pin_file:
                with tempfile.NamedTemporaryFile(delete=False, suffix="pin.txt") as tmp:
                    tmp.write(self.content_pin)
                    pin_file = tmp.name
                    tmp_pin_file = tmp.name
        else:
            env_file = env_file.get_path_or_uri()
            deploy_file = self.post_deploy_file
            pin_file = self.pin_file

        env_path = self.address

        if self.is_containerized:
            if not dryrun:
                try:
                    shell.check_output(
                        singularity.shellcmd(
                            self._container_img.path,
                            "[ -d '{}' ]".format(env_path),
                            args=self._singularity_args,
                            envvars=self.get_singularity_envvars(),
                            quiet=True,
                        ),
                        stderr=subprocess.PIPE,
                    )
                except subprocess.CalledProcessError as e:
                    raise WorkflowError(
                        "Unable to find environment in container image. "
                        "Maybe a conda environment was modified without containerizing again "
                        "(see snakemake --containerize)?\nDetails:\n{}\n{}".format(
                            e, e.stderr.decode()
                        )
                    )
                return env_path
            else:
                # env should be present in the container
                return env_path

        # Check for broken environment
        if os.path.exists(
            os.path.join(env_path, "env_setup_start")
        ) and not os.path.exists(os.path.join(env_path, "env_setup_done")):
            if dryrun:
                logger.info(
                    "Incomplete Conda environment {} will be recreated.".format(
                        self.file.simplify_path()
                    )
                )
            else:
                logger.info(
                    "Removing incomplete Conda environment {}...".format(
                        self.file.simplify_path()
                    )
                )
                shutil.rmtree(env_path, ignore_errors=True)

        # Create environment if not already present.
        if not os.path.exists(env_path):
            if dryrun:
                logger.info(
                    "Conda environment {} will be created.".format(
                        self.file.simplify_path()
                    )
                )
                return env_path
            logger.info(
                "Creating conda environment {}...".format(self.file.simplify_path())
            )
            env_archive = self.archive_file
            try:
                # Touch "start" flag file
                os.makedirs(env_path, exist_ok=True)
                with open(os.path.join(env_path, "env_setup_start"), "a") as f:
                    pass

                # Check if env archive exists. Use that if present.
                if os.path.exists(env_archive):
                    logger.info("Installing archived conda packages.")
                    pkg_list = os.path.join(env_archive, "packages.txt")
                    if os.path.exists(pkg_list):
                        # read pacakges in correct order
                        # this is for newer env archives where the package list
                        # was stored
                        packages = [
                            os.path.join(env_archive, pkg.rstrip())
                            for pkg in open(pkg_list)
                        ]
                    else:
                        # guess order
                        packages = glob(os.path.join(env_archive, "*.tar.bz2"))

                    # install packages manually from env archive
                    cmd = " ".join(
                        [
                            "conda",
                            "create",
                            "--quiet",
                            "--yes",
                            "--prefix '{}'".format(env_path),
                        ]
                        + packages
                    )
                    if self._container_img:
                        cmd = singularity.shellcmd(
                            self._container_img.path,
                            cmd,
                            args=self._singularity_args,
                            envvars=self.get_singularity_envvars(),
                        )
                    out = shell.check_output(
                        cmd, stderr=subprocess.STDOUT, universal_newlines=True
                    )
                else:

                    def create_env(env_file, filetype="yaml"):
                        # Copy env file to env_path (because they can be on
                        # different volumes and singularity should only mount one).
                        # In addition, this allows to immediately see what an
                        # environment in .snakemake/conda contains.
                        target_env_file = env_path + f".{filetype}"
                        shutil.copy(env_file, target_env_file)

                        logger.info("Downloading and installing remote packages.")

                        strict_priority = (
                            ["conda config --set channel_priority strict &&"]
                            if self._container_img
                            else []
                        )

                        subcommand = [self.frontend]
                        yes_flag = ["--yes"]
                        if filetype == "yaml":
                            subcommand.append("env")
                            yes_flag = []

                        cmd = (
                            strict_priority
                            + subcommand
                            + [
                                "create",
                                "--quiet",
                                '--file "{}"'.format(target_env_file),
                                '--prefix "{}"'.format(env_path),
                            ]
                            + yes_flag
                        )
                        cmd = " ".join(cmd)
                        if self._container_img:
                            cmd = singularity.shellcmd(
                                self._container_img.path,
                                cmd,
                                args=self._singularity_args,
                                envvars=self.get_singularity_envvars(),
                            )
                        out = shell.check_output(
                            cmd, stderr=subprocess.STDOUT, universal_newlines=True
                        )

                        # cleanup if requested
                        if self._cleanup is CondaCleanupMode.tarballs:
                            logger.info("Cleaning up conda package tarballs.")
                            shell.check_output("conda clean -y --tarballs")
                        elif self._cleanup is CondaCleanupMode.cache:
                            logger.info(
                                "Cleaning up conda package tarballs and package cache."
                            )
                            shell.check_output("conda clean -y --tarballs --packages")
                        return out

                    if pin_file is not None:
                        try:
                            logger.info(
                                f"Using pinnings from {self.pin_file.get_path_or_uri()}."
                            )
                            out = create_env(pin_file, filetype="pin.txt")
                        except subprocess.CalledProcessError as e:
                            # remove potential partially installed environment
                            shutil.rmtree(env_path, ignore_errors=True)
                            advice = ""
                            if isinstance(self.file, LocalSourceFile):
                                advice = (
                                    " If that works, make sure to update the pin file with "
                                    f"'snakedeploy pin-conda-env {self.file.get_path_or_uri()}'."
                                )
                            logger.warning(
                                f"Failed to install conda environment from pin file ({self.pin_file.get_path_or_uri()}). "
                                f"Trying regular environment definition file.{advice}"
                            )
                            out = create_env(env_file, filetype="yaml")
                    else:
                        out = create_env(env_file, filetype="yaml")

                # Execute post-deplay script if present
                if deploy_file:
                    target_deploy_file = env_path + ".post-deploy.sh"
                    shutil.copy(deploy_file, target_deploy_file)
                    self.execute_deployment_script(env_file, target_deploy_file)

                # Touch "done" flag file
                with open(os.path.join(env_path, "env_setup_done"), "a") as f:
                    pass

                logger.debug(out)
                logger.info(
                    f"Environment for {self.file.get_path_or_uri()} created (location: {os.path.relpath(env_path)})"
                )
            except subprocess.CalledProcessError as e:
                # remove potential partially installed environment
                shutil.rmtree(env_path, ignore_errors=True)
                raise CreateCondaEnvironmentException(
                    f"Could not create conda environment from {env_file}:\nCommand:\n{e.cmd}\nOutput:\n{e.output}"
                )

        if tmp_env_file:
            # temporary file was created
            os.remove(tmp_env_file)
        if tmp_deploy_file:
            os.remove(tmp_deploy_file)

        return env_path

    @classmethod
    def get_singularity_envvars(self):
        return {"CONDA_PKGS_DIRS": "/tmp/conda/{}".format(uuid.uuid4())}

    def __hash__(self):
        # this hash is only for object comparison, not for env paths
        if self.is_named:
            return hash(self.name)
        else:
            return hash(self.file)

    def __eq__(self, other):
        if isinstance(other, Env):
            if self.is_named:
                return self.name == other.name
            else:
                return self.file == other.file
        return False


class Conda:
    instances = dict()
    lock = threading.Lock()

    def __new__(cls, container_img=None, prefix_path=None):
        with cls.lock:
            if container_img not in cls.instances:
                inst = super().__new__(cls)
                cls.instances[container_img] = inst
                return inst
            else:
                return cls.instances[container_img]

    def __init__(self, container_img=None, prefix_path=None):
        if not self.is_initialized:  # avoid superfluous init calls
            from snakemake.deployment import singularity
            from snakemake.shell import shell

            if isinstance(container_img, singularity.Image):
                container_img = container_img.path
            self.container_img = container_img

            self.info = json.loads(
                shell.check_output(
                    self._get_cmd("conda info --json"), universal_newlines=True
                )
            )

            if prefix_path is None or container_img is not None:
                self.prefix_path = self.info["conda_prefix"]
            else:
                self.prefix_path = prefix_path

            self.platform = self.info["platform"]

            # check conda installation
            self._check()

    @property
    def is_initialized(self):
        return hasattr(self, "prefix_path")

    def _get_cmd(self, cmd):
        if self.container_img:
            return singularity.shellcmd(self.container_img, cmd, quiet=True)
        return cmd

    def _check(self):
        from snakemake.shell import shell

        # Use type here since conda now is a function.
        # type allows to check for both functions and regular commands.
        if not ON_WINDOWS or shell.get_executable():
            locate_cmd = "type conda"
        else:
            locate_cmd = "where conda"

        try:
            shell.check_output(self._get_cmd(locate_cmd), stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            if self.container_img:
                raise CreateCondaEnvironmentException(
                    "The 'conda' command is not "
                    "available inside "
                    "your singularity container "
                    "image. Snakemake mounts "
                    "your conda installation "
                    "into singularity. "
                    "Sometimes, this can fail "
                    "because of shell restrictions. "
                    "It has been tested to work "
                    "with docker://ubuntu, but "
                    "it e.g. fails with "
                    "docker://bash "
                )
            else:
                raise CreateCondaEnvironmentException(
                    "The 'conda' command is not "
                    "available in the "
                    "shell {} that will be "
                    "used by Snakemake. You have "
                    "to ensure that it is in your "
                    "PATH, e.g., first activating "
                    "the conda base environment "
                    "with `conda activate base`.".format(shell.get_executable())
                )
        try:
            self._check_version()
            self._check_condarc()
        except subprocess.CalledProcessError as e:
            raise CreateCondaEnvironmentException(
                "Unable to check conda installation:\n" + e.stderr.decode()
            )

    def _check_version(self):
        from snakemake.shell import shell

        version = shell.check_output(
            self._get_cmd("conda --version"),
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        version_matches = re.findall("\d+.\d+.\d+", version)
        if len(version_matches) != 1:
            raise WorkflowError(
                f"Unable to determine conda version. 'conda --version' returned {version}"
            )
        else:
            version = version_matches[0]
        if StrictVersion(version) < StrictVersion("4.2"):
            raise CreateCondaEnvironmentException(
                "Conda must be version 4.2 or later, found version {}.".format(version)
            )

    def _check_condarc(self):
        if self.container_img:
            # Do not check for strict priorities when running conda in an image
            # Instead, we set priorities to strict ourselves in the image.
            return
        from snakemake.shell import shell

        res = json.loads(
            shell.check_output(
                self._get_cmd("conda config --get channel_priority --json"),
                universal_newlines=True,
                stderr=subprocess.PIPE,
            )
        )
        if res["get"].get("channel_priority") != "strict":
            raise CreateCondaEnvironmentException(
                "Your conda installation is not configured to use strict channel priorities. "
                "This is however crucial for having robust and correct environments (for details, "
                "see https://conda-forge.org/docs/user/tipsandtricks.html). "
                "Please configure strict priorities by executing 'conda config --set channel_priority strict'."
            )

    def bin_path(self):
        if ON_WINDOWS:
            return os.path.join(self.prefix_path, "Scripts")
        else:
            return os.path.join(self.prefix_path, "bin")

    def shellcmd(self, env_address, cmd):
        # get path to activate script
        activate = os.path.join(self.bin_path(), "activate")

        if ON_WINDOWS:
            activate = activate.replace("\\", "/")
            env_address = env_address.replace("\\", "/")

        return "source {} '{}'; {}".format(activate, env_address, cmd)

    def shellcmd_win(self, env_address, cmd):
        """Prepend the windows activate bat script."""
        # get path to activate script
        activate = os.path.join(self.bin_path(), "activate.bat").replace("\\", "/")
        env_address = env_address.replace("\\", "/")

        return '"{}" "{}"&&{}'.format(activate, env_address, cmd)


def is_mamba_available():
    return shutil.which("mamba") is not None


class CondaEnvSpec(ABC):
    @abstractmethod
    def apply_wildcards(self, wildcards):
        ...

    @abstractmethod
    def get_conda_env(self, workflow, env_dir=None, container_img=None, cleanup=None):
        ...

    @abstractmethod
    def check(self):
        ...

    @property
    def is_file(self):
        return False

    @property
    @abstractmethod
    def contains_wildcard(self):
        ...

    @abstractmethod
    def __hash__(self):
        ...

    @abstractmethod
    def __eq__(self, other):
        ...


class CondaEnvFileSpec(CondaEnvSpec):
    def __init__(self, filepath, rule=None):
        if isinstance(filepath, SourceFile):
            self.file = IOFile(str(filepath.get_path_or_uri()), rule=rule)
        elif isinstance(filepath, _IOFile):
            self.file = filepath
        else:
            self.file = IOFile(filepath, rule=rule)

    def apply_wildcards(self, wildcards, rule):
        filepath = self.file.apply_wildcards(wildcards)
        if is_local_file(filepath):
            # Normalize 'file:///my/path.yml' to '/my/path.yml'
            filepath = parse_uri(filepath).uri_path
        return CondaEnvFileSpec(filepath, rule)

    def check(self):
        self.file.check()

    def get_conda_env(self, workflow, env_dir=None, container_img=None, cleanup=None):
        return Env(
            workflow,
            env_file=self.file,
            env_dir=env_dir,
            container_img=container_img,
            cleanup=cleanup,
        )

    @property
    def is_file(self):
        return True

    @property
    def contains_wildcard(self):
        return contains_wildcard(self.file)

    def __hash__(self):
        return hash(self.file)

    def __eq__(self, other):
        return self.file == other.file


class CondaEnvNameSpec(CondaEnvSpec):
    def __init__(self, name: str):
        self.name = name

    def apply_wildcards(self, wildcards, _):
        return CondaEnvNameSpec(apply_wildcards(self.name, wildcards))

    def get_conda_env(self, workflow, env_dir=None, container_img=None, cleanup=None):
        return Env(
            workflow,
            env_name=self.name,
            env_dir=env_dir,
            container_img=container_img,
            cleanup=cleanup,
        )

    def check(self):
        # not a file, nothing to check here
        pass

    @property
    def contains_wildcard(self):
        return contains_wildcard(self.name)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name


def is_conda_env_file(spec):
    if isinstance(spec, SourceFile):
        spec = spec.get_filename()

    return spec.endswith(".yaml") or spec.endswith(".yml")
