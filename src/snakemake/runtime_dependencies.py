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

from snakemake import __version__


@dataclass(frozen=True)
class Package:
    """Represents a Python package with an optional version."""

    name: str
    version: Optional[str]

    @classmethod
    def from_str(cls, pkg_str: str) -> Self:
        if "==" in pkg_str:
            name, version = pkg_str.split("==", 1)
            return cls(name, version)
        else:
            return cls(pkg_str, None)

    def __str__(self) -> str:
        return f"{self.name}=={self.version}" if self.version else self.name


class PackageType(Enum):
    """Represents the type of a runtime dependency package."""

    PLUGIN = 0
    WORKFLOW = 1


class RuntimeDependencyManager:
    """Manages workflow runtime dependencies for Snakemake, including plugin and
    auxiliary packages.
    """

    _instances: Self = None
    _lock_lifetime = timedelta(seconds=30)

    def __new__(cls) -> Self:
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.__init__()
        return cls._instance

    def __init__(self):
        if hasattr(self, "_packages"):
            # instantiated before, skip
            return
        self._prefixes: Map[PackageType, Path] = {
            PackageType.PLUGIN: platformdirs.user_cache_path(
                appname="snakemake", version=__version__, ensure_exists=True
            )
            / "plugins",
        }
        self._packages: Dict[PackageType, Set[Package]] = {
            PackageType.WORKFLOW: set(),
            PackageType.PLUGIN: set(),
        }

    def update_workflow_prefix(self) -> None:
        self._prefixes[PackageType.WORKFLOW] = Path.cwd() / ".snakemake" / "workflow_dependencies"

    def add_plugin_package(self, name: str) -> None:
        self._add_package(PackageType.PLUGIN, name)

    def add_workflow_package(self, name: str, version: Optional[str] = None) -> None:
        self._add_package(
            PackageType.WORKFLOW,
            name,
            version,
        )

    def deploy_plugin_packages(self) -> None:
        """Deploys the runtime dependencies for all package types, ensuring that
        plugin packages are deployed before auxiliary packages.
        """
        if not self._packages[PackageType.PLUGIN]:
            return
        # retrieve packages from current environment (e.g. snakemake)
        prior_packages = set(get_packages_in_current_env())
        # deploy plugin packages, considering the packages from the current
        # environment as prior packages to ensure compatibility
        self._deploy_packages_per_type(PackageType.PLUGIN, prior_packages)

    def deploy_workflow_packages(self) -> None:
        if not self._packages[PackageType.WORKFLOW]:
            return

        # retrieve packages from current environment (e.g. snakemake)
        prior_env_packages = set(get_packages_in_current_env())

        # get all packages and versions installed under the plugin package prefix
        # (including other packages installed there before)
        prior_plugin_packages = set(get_packages_in_prefix(self._prefixes[PackageType.PLUGIN]))

        # deploy auxiliary packages, considering the plugin packages as additional
        # prior packages to ensure compatibility
        self._deploy_packages_per_type(PackageType.WORKFLOW, prior_env_packages, prior_plugin_packages)

    def _add_package(
        self, package_type: PackageType, name: str, version: Optional[str] = None
    ) -> None:
        """Adds a package to the set of runtime dependencies for the given
        package type.
        """
        package = Package(name, version)
        self._packages[package_type].add(package)

    def _deploy_packages_per_type(
        self, package_type: PackageType, *prior_package_sets: Set[Package]
    ) -> None:
        """Deploys the runtime dependencies for the given package type, optionally
        considering additional packages (e.g. plugin packages when deploying
        auxiliary packages).
        """
        prior_packages = {}
        for pkg_set in prior_package_sets:
            for pkg in pkg_set:
                prior_packages[pkg.name] = pkg

        prefix = self._prefixes[package_type]
        lock = self._lock(prefix)
        with lock:
            try:
                refresher_thread = Thread(
                    target=partial(self._refresh_lock, lock), daemon=True
                )
                refresher_thread.start()
                # add the actually requested packages (may overwrite versions)
                for pkg in self._packages[package_type]:
                    prior_packages[pkg.name] = pkg

                # let uv solve the requested packages together with the additional
                # packages and the one from the current environment,
                # thereby ensuring that the resulting set of packages in consistent
                # and compatible versions are installed.
                res = sp.run(
                    ["uv", "pip", "install", "--prefix", prefix, "--dry-run"]
                    + prior_packages,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                # extract the requested packages with their determined versions from the
                # uv output
                posterior_packages = set(
                    package
                    for package in parse_uv_pip_dry_run_output(res.stderr)
                    if package.name in self._packages[package_type]
                )
                # explicitly install the determined versions of the requested
                # packages under the prefix
                sp.run(
                    ["uv", "pip", "install", "--prefix", self.prefix]
                    + [str(pkg) for pkg in posterior_packages],
                    check=True,
                )
            finally:
                refresher_thread.join()

    def _lock(self, prefix: Path) -> Lock:
        lockfile = prefix / "deploy.lock"
        lock = Lock(lockfile)
        lock.lifetime = self._lock_lifetime
        return lock

    @classmethod
    def _refresh_lock(cls, lock: Lock) -> None:
        while lock.is_locked:
            lock.refresh()
            time.sleep(cls._lock_lifetime.total_seconds() / 2)


def parse_uv_pip_dry_run_output(output: str) -> Iterable[Package]:
    """Parses the output of `uv pip install --dry-run` to extract the packages that would be installed."""
    install_prefix = " + "
    for line in output.splitlines():
        line = line.strip()
        if line.startswith(install_prefix):
            pkg_str = line.removeprefix(install_prefix)
            yield Package.from_str(pkg_str)


def get_packages_in_prefix(prefix: Path) -> Iterable[Package]:
    """Retrieves the packages installed under the given prefix using
    `uv pip list --prefix <prefix>`.
    """
    return map(
        Package,
        json.loads(
            sp.run(
                ["uv", "pip", "list", "--prefix", prefix],
                check=True,
                capture_output=True,
            ).stdout
        ),
    )


def get_packages_in_current_env() -> Iterable[Package]:
    """Retrieves the packages installed in the current environment
    using `uv pip list`.
    """
    python_exec = Path(shutil.which("python"))
    prefix = python_exec.parent.parent
    return get_packages_in_prefix(prefix)
