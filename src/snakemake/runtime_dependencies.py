from itertools import chain
import threading
import sys
from typing import Dict
from immutables import Map
from snakemake_interface_common.exceptions import WorkflowError
import shutil
import json
from enum import Enum
import platformdirs
from functools import partial
from typing import Iterable
from typing import Self
from threading import Thread
from datetime import timedelta
from pathlib import Path
from typing import Set
from dataclasses import dataclass
from typing import Optional
import subprocess as sp
import time

from flufl.lock import Lock
from packaging.requirements import Requirement

from snakemake import __version__


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
        self._prefixes: Map[PackageType, Path] = {
            PackageType.GLOBAL: platformdirs.user_cache_path(
                appname="snakemake", version=__version__, ensure_exists=True
            )
            / "global_dependencies",
            PackageType.WORKFLOW: deployment_prefix / "workflow_dependencies",
        }
        self._packages: Dict[PackageType, Dict[str, Requirement]] = {
            PackageType.WORKFLOW: dict(),
            PackageType.GLOBAL: dict(),
        }

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
        prior_env_packages = set(get_packages_in_current_env())

        try:
            # deploy global packages, considering the packages from the current
            # environment as prior packages to ensure compatibility
            self._deploy_packages_per_type(
                self._packages[PackageType.GLOBAL].values(),
                self._prefixes[PackageType.GLOBAL],
                prior_env_packages,
            )
        except sp.CalledProcessError as e:
            raise WorkflowError(
                f"Failed to deploy auxilliary global python packages (--with-pkgs option): {e.stderr}"
            )

        # get final versions of global packages and all their dependencies in the prefix
        prior_global_packages = set(
            get_packages_in_prefix(self._prefixes[PackageType.GLOBAL])
        )

        try:
            self._deploy_packages_per_type(
                self._packages[PackageType.GLOBAL].values(),
                self._prefixes[PackageType.WORKFLOW],
                prior_env_packages,
                prior_global_packages,
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
                )
            except sp.CalledProcessError as e2:
                raise WorkflowError(
                    "Failed to deploy auxilliary workflow python packages (--workflow-with-pkgs option)."
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
    ) -> None:
        """Deploys the runtime dependencies for the given package type, optionally
        considering additional packages (e.g. plugin packages when deploying
        auxiliary packages).
        """
        packages = list(packages)
        if not packages:
            # no packages requested, stop early
            return

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
                    ["uv", "pip", "install", "--target", str(prefix), "--dry-run"]
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
                )
                if posterior_packages:
                    # explicitly install the determined versions of the requested
                    # packages under the prefix
                    sp.run(
                        ["uv", "pip", "install", "--no-deps", "--target", str(prefix)]
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


def parse_uv_pip_dry_run_output(output: str) -> Iterable[Requirement]:
    """Parses the output of `uv pip install --dry-run` to extract the packages that would be installed."""
    install_prefix = "+ "
    for line in output.splitlines():
        line = line.strip()
        if line.startswith(install_prefix):
            pkg_str = line.removeprefix(install_prefix)
            yield Requirement(pkg_str)


def get_packages_in_prefix(prefix: Path) -> Iterable[Requirement]:
    """Retrieves the packages installed under the given prefix using
    `uv pip list --prefix <prefix>`.
    """

    def entry_to_req(entry: Dict[str, str]) -> Requirement:
        extras = "'".join(entry.get("extras", []))
        if extras:
            extras = f"[{extras}]"
        return Requirement(f"{entry['name']}{extras}=={entry['version']}")

    return map(
        entry_to_req,
        json.loads(
            sp.run(
                ["uv", "pip", "list", "--prefix", prefix, "--format", "json"],
                check=True,
                capture_output=True,
            ).stdout
        ),
    )


def get_packages_in_current_env() -> Iterable[Requirement]:
    """Retrieves the packages installed in the current environment
    using `uv pip list`.
    """
    python_exec = Path(sys.executable)
    prefix = python_exec.parent.parent
    return get_packages_in_prefix(prefix)
