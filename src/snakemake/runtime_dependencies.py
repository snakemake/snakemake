import hashlib
from itertools import chain
import threading
import sys
from typing import Dict

from importlib.metadata import Distribution
from immutables import Map
from snakemake.common.constants import ON_WINDOWS
from snakemake_interface_common.exceptions import WorkflowError
import shutil
import json
from enum import Enum
import platformdirs
from functools import partial
from typing import Iterable
from threading import Thread
from datetime import timedelta
from pathlib import Path
from typing import Set
from typing import Optional
import subprocess as sp
import time

from flufl.lock import Lock
from packaging.requirements import Requirement

from snakemake import __version__
from snakemake.logging import logger
from snakemake.common.misc import get_minor_version


class PackageType(Enum):
    """Represents the type of a runtime dependency package."""

    GLOBAL = 0
    WORKFLOW = 1


class RuntimeDependencyManager:
    """Manages workflow runtime dependencies for Snakemake, including plugin and
    auxiliary packages.
    """

    _lock_lifetime = timedelta(seconds=30)

    def __init__(self, deployment_prefix: Path):
        version_path = f"py_{sys.version_info.major}.{sys.version_info.minor}/snakemake_{get_minor_version()}"
        if is_editable():
            version_path += "+edits"
        self._prefixes: Map[PackageType, Path] = {
            PackageType.GLOBAL: platformdirs.user_cache_path(
                appname="snakemake", version=version_path, ensure_exists=True
            )
            / "global_dependencies",
            PackageType.WORKFLOW: deployment_prefix
            / "workflow_dependencies"
            / version_path,
        }
        self._packages: Dict[PackageType, Dict[str, Requirement]] = {
            PackageType.WORKFLOW: dict(),
            PackageType.GLOBAL: dict(),
        }
        self._uv_exec: Optional[str] = None

    @property
    def uv_exec(self) -> str:
        if self._uv_exec is not None:
            return self._uv_exec

        suffix = ".exe" if ON_WINDOWS else ""
        if (in_env_uv := Path(sys.executable).parent / f"uv{suffix}").exists():
            self._uv_exec = str(in_env_uv)
        elif (global_uv := shutil.which("uv")) is not None:
            self._uv_exec = global_uv
        else:
            raise WorkflowError("Unable to find 'uv' executable in the environment.")

        return self._uv_exec

    def update_workflow_prefix(self) -> None:
        self._prefixes[PackageType.WORKFLOW] = (
            Path.cwd() / ".snakemake" / "workflow_dependencies"
        )

    def add_global_packages(self, *pkgs: str) -> None:
        for pkg in pkgs:
            self._add_package(PackageType.GLOBAL, Requirement(pkg))

    def add_workflow_packages(self, *pkgs: str) -> None:
        for pkg in pkgs:
            self._add_package(
                PackageType.WORKFLOW,
                Requirement(pkg),
            )

    def deploy_packages(self) -> None:
        # retrieve packages from current environment (e.g. snakemake)
        prior_env_packages = set(self.get_packages_in_current_env())

        try:
            # deploy global packages, considering the packages from the current
            # environment as prior packages to ensure compatibility
            self._deploy_packages_per_type(
                self._packages[PackageType.GLOBAL].values(),
                self._prefixes[PackageType.GLOBAL],
                prior_env_packages,
                separate_by_hash=False,
            )
        except sp.CalledProcessError as e:
            raise WorkflowError(
                f"Failed to deploy auxiliary global python packages (--with-pkgs option): {e.stderr}"
            )

        # get final versions of global packages and all their dependencies in the prefix
        prior_global_packages = set(
            self.get_packages_in_prefix(self._prefixes[PackageType.GLOBAL])
        )

        try:
            self._deploy_packages_per_type(
                self._packages[PackageType.WORKFLOW].values(),
                self._prefixes[PackageType.WORKFLOW],
                prior_env_packages,
                prior_global_packages,
                # workflow specific packages have to be separated by their hash to avoid conflicts
                # between workflows and between different profiles of the same workflow
                separate_by_hash=True,
            )
        except sp.CalledProcessError as e:
            # try to solve together as a fallback
            try:
                self._deploy_packages_per_type(
                    chain(
                        self._packages[PackageType.GLOBAL].values(),
                        self._packages[PackageType.WORKFLOW].values(),
                    ),
                    self._prefixes[PackageType.WORKFLOW],
                    prior_env_packages,
                    separate_by_hash=True,
                )
            except sp.CalledProcessError as e2:
                raise WorkflowError(
                    "Failed to deploy auxiliary workflow python packages (--workflow-with-pkgs option)."
                    f"\nSeparate solve: {e.stderr}\nJoint fallback solve: {e2.stderr}"
                )

    def _add_package(self, package_type: PackageType, pkg: Requirement) -> None:
        """Adds a package to the set of runtime dependencies for the given
        package type.
        """
        self._packages[package_type][pkg.name] = pkg

    def _deploy_packages_per_type(
        self,
        packages: Iterable[Requirement],
        prefix: Path,
        *prior_package_sets: Set[Requirement],
        separate_by_hash: bool = False,
    ) -> None:
        """Deploys the runtime dependencies for the given package type, optionally
        considering additional packages (e.g. plugin packages when deploying
        auxiliary packages).
        """
        packages = list(packages)
        if not packages:
            # no packages requested, stop early
            return

        if separate_by_hash:
            prefix = prefix / hash_packages(packages)

        prior_packages = {}
        for pkg_set in prior_package_sets:
            for pkg in pkg_set:
                prior_packages[pkg.name] = pkg

        # add the actually requested packages (may overwrite versions)
        for pkg in packages:
            prior_packages[pkg.name] = pkg

        requested_pkg_names = {pkg.name for pkg in packages}

        lock = self._lock(prefix)
        with lock:
            stop_refresher = threading.Event()
            refresher_thread = Thread(
                target=partial(self._refresh_lock, lock, stop_refresher), daemon=True
            )
            refresher_thread.start()
            try:
                # let uv solve the requested packages together with the additional
                # packages and the one from the current environment,
                # thereby ensuring that the resulting set of packages in consistent
                # and compatible versions are installed.
                res = sp.run(
                    [
                        self.uv_exec,
                        "pip",
                        "install",
                        "--target",
                        str(prefix),
                        "--dry-run",
                        "--python",
                        sys.executable,
                    ]
                    + [str(pkg) for pkg in prior_packages.values()],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                # extract the requested packages with their determined versions from the
                # uv output
                posterior_packages = set(
                    package
                    for package in parse_uv_pip_dry_run_output(res.stderr)
                    if package.name in requested_pkg_names
                    or package.name not in prior_packages
                )
                if posterior_packages:
                    print(
                        f"Deploying packages (prefix: {prefix}): {', '.join(str(pkg) for pkg in posterior_packages)}",
                        file=sys.stderr,
                    )
                    # remove old packages that are not needed anymore.
                    # get all packages in the prefix and remove those that are not in the dryrun result
                    to_uninstall = [
                        str(pkg)
                        for pkg in self.get_packages_in_prefix(prefix)
                        if pkg not in posterior_packages
                    ]
                    if to_uninstall:
                        sp.run(
                            [
                                self.uv_exec,
                                "pip",
                                "uninstall",
                                "-y",
                                "--target",
                                str(prefix),
                                "--python",
                                sys.executable,
                            ]
                            + to_uninstall,
                            check=True,
                            capture_output=True,
                            text=True,
                        )
                    # install the new packages that are needed according to the dryrun result
                    sp.run(
                        [
                            self.uv_exec,
                            "pip",
                            "install",
                            "--no-deps",  # we have an explicitl list, the rest is in the prior packages
                            "--target",
                            str(prefix),
                            "--python",
                            sys.executable,
                        ]
                        + [str(pkg) for pkg in posterior_packages],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                sys.path.insert(0, str(prefix))
            finally:
                stop_refresher.set()

    def _lock(self, prefix: Path) -> Lock:
        prefix.mkdir(parents=True, exist_ok=True)
        lockfile = prefix / "deploy.lock"
        lock = Lock(str(lockfile))
        lock.lifetime = self._lock_lifetime
        return lock

    @classmethod
    def _refresh_lock(cls, lock: Lock, stop: threading.Event) -> None:
        while not stop.is_set():
            lock.refresh()
            time.sleep(cls._lock_lifetime.total_seconds() / 2)

    def get_packages_in_prefix(self, prefix: Path) -> Iterable[Requirement]:
        """Retrieves the packages installed under the given prefix using
        `uv pip list --prefix <prefix>`.
        """

        def entry_to_req(entry: Dict[str, str]) -> Requirement:
            extras = "'".join(entry.get("extras", []))
            if extras:
                extras = f"[{extras}]"
            if editable_location := entry.get("editable_project_location"):
                return Requirement(
                    f"{entry['name']}{extras} @ {Path(editable_location).as_uri()}"
                )
            return Requirement(f"{entry['name']}{extras}=={entry['version']}")

        return map(
            entry_to_req,
            json.loads(
                sp.run(
                    [
                        self.uv_exec,
                        "pip",
                        "list",
                        "--prefix",
                        prefix,
                        "--format",
                        "json",
                    ],
                    check=True,
                    capture_output=True,
                ).stdout
            ),
        )

    def get_packages_in_current_env(self) -> Iterable[Requirement]:
        """Retrieves the packages installed in the current environment
        using `uv pip list`.
        """
        python_exec = Path(sys.executable)
        prefix = python_exec.parent.parent
        return self.get_packages_in_prefix(prefix)


def parse_uv_pip_dry_run_output(output: str) -> Iterable[Requirement]:
    """Parses the output of `uv pip install --dry-run` to extract the packages that would be installed."""
    install_prefix = "+ "
    for line in output.splitlines():
        line = line.strip()
        if line.startswith(install_prefix):
            pkg_str = line.removeprefix(install_prefix)
            yield Requirement(pkg_str)


def is_editable() -> bool:
    direct_url = Distribution.from_name("snakemake").read_text("direct_url.json")
    if direct_url is None:
        return False
    dir_info = json.loads(direct_url).get("dir_info", {})
    return dir_info.get("editable", False)


def hash_packages(packages: Iterable[Requirement]) -> str:
    """Generates a hash for the given set of packages."""
    m = hashlib.sha256()
    for pkg in sorted(packages, key=str):
        m.update(str(pkg).encode())
    return m.hexdigest()
