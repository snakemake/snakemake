from snakemake.logging import logger
from snakemake.settings.enums import Quietness
from snakemake_interface_executor_plugins.settings import SharedFSUsage
from snakemake.common import is_local_file
from copy import copy
from typing import List
from typing import Callable
from typing import Union
from typing import Optional
from dataclasses import dataclass
from snakemake.common import get_snakemake_searchpaths
from pathlib import Path
import tempfile
from typing import Any, Dict, Iterable
import os

import snakemake
from snakemake_interface_common.exceptions import WorkflowError
from snakemake_interface_software_deployment_plugins.registry import (
    SoftwareDeploymentPluginRegistry,
)
from snakemake_interface_software_deployment_plugins import (
    EnvSpecBase,
    EnvBase,
    EnvSpecSourceFile,
    ShellExecutable,
)
from snakemake.common import get_function_params, overwrite_function_params
from snakemake_software_deployment_plugin_conda import EnvSpec as CondaEnvSpec
from snakemake_software_deployment_plugin_container import EnvSpec as ContainerEnvSpec
from snakemake_software_deployment_plugin_envmodules import EnvSpec as EnvModuleEnvSpec


class SoftwareDeploymentManager:
    def __init__(self, workflow: "snakemake.workflow.Workflow"):
        from snakemake.shell import shell
        from snakemake.workflow import Workflow

        self.workflow: Workflow = workflow
        self.specs_to_envs = {}
        self.env_instances = {}
        self.shell_executable = ShellExecutable(
            executable=shell.get_executable() or "/bin/sh",
            command_arg="-c",
        )
        self.registry = SoftwareDeploymentPluginRegistry()
        self.plugins = {}
        for plugin_name, plugin in self.registry.plugins.items():
            if plugin_name not in self.workflow.deployment_settings.deployment_methods:
                continue
            kind = plugin.common_settings.provides
            if kind not in self.plugins:
                self.plugins[kind] = plugin
            else:
                raise WorkflowError(
                    "Multiple plugins for the same kind are not allowed. "
                    f"Choose one of {plugin_name} and {self.plugins[kind].name}."
                )

    def collect_envs(self, jobs: Iterable["snakemake.jobs.Job"]) -> None:
        shared_envs = (
            SharedFSUsage.SOFTWARE_DEPLOYMENT
            in self.workflow.storage_settings.shared_fs_usage
            or SharedFSUsage.SOFTWARE_DEPLOYMENT_CACHE
            in self.workflow.storage_settings.shared_fs_usage
        )
        for job in jobs:
            if (
                job.software_env_spec is not None
                and job.software_env_spec not in self.specs_to_envs
                and (
                    job.is_local
                    or shared_envs
                    or (
                        self.workflow.remote_exec
                        and SharedFSUsage.SOFTWARE_DEPLOYMENT
                        not in self.workflow.storage_settings.shared_fs_usage
                    )
                )
            ):
                self.get_env(job.software_env_spec)

    async def cache_envs(self, jobs: Iterable["snakemake.jobs.Job"]) -> None:
        shared_cache = (
            SharedFSUsage.SOFTWARE_DEPLOYMENT_CACHE
            in self.workflow.storage_settings.shared_fs_usage
        )
        if not shared_cache and not self.workflow.remote_exec:
            return

        assets = {}
        for job in jobs:
            env = job.software_env
            if env is not None and env.is_cacheable():
                for asset in await env.get_cache_assets():
                    if (env.cache_path / asset).exists():
                        continue
                    # if multiple envs share the same asset, we can simply cache
                    # it for one of them, so we do not check for duplicates here.
                    assets[asset] = env

        if self.workflow.dryrun and assets:
            logger.info(
                f"Will cache {len(assets)} software environment assets.",
                extra={"quietness": Quietness.RULES},
            )
        elif assets:
            logger.info(
                f"Caching {len(assets)} software environment assets.",
                extra={"quietness": Quietness.RULES},
            )

            for asset, env in assets.items():
                await env.managed_cache_asset(asset)

    async def deploy_envs(self, jobs: Iterable["snakemake.jobs.Job"]) -> None:
        shared_deploy = (
            SharedFSUsage.SOFTWARE_DEPLOYMENT
            in self.workflow.storage_settings.shared_fs_usage
        )
        if not shared_deploy and not self.workflow.remote_exec:
            return

        envs = set()
        for job in jobs:
            env = job.software_env
            if (
                env is not None
                and env.is_deployable()
                and not env.deployment_path.exists()
            ):
                envs.add(env)

        if self.workflow.dryrun and envs:
            logger.info(
                f"Will deploy {len(envs)} software environments.",
                extra={"quietness": Quietness.RULES},
            )
        elif envs:
            logger.info(
                f"Deploying {len(envs)} software environments.",
                extra={"quietness": Quietness.RULES},
            )
            for env in envs:
                await env.managed_deploy()

    async def cleanup_envs(self, jobs: Iterable["snakemake.jobs.Job"]) -> None:
        for job in jobs:
            env = job.software_env
            if env is not None:
                await env.cleanup_cache()
                env.managed_remove()

    def get_env(
        self, env_spec: EnvSpecBase, mountpoints: List[Path] = []
    ) -> Optional[EnvBase]:
        if env_spec in self.specs_to_envs:
            return self.specs_to_envs[env_spec]

        if env_spec.kind not in self.plugins:
            if env_spec.fallback is not None:
                return self.get_env(env_spec.fallback)
            else:
                # no method activated that can yield an env here
                return None

        deployment_prefix = (
            self.workflow.deployment_settings.deployment_prefix / env_spec.kind
        )
        cache_prefix = self.workflow.deployment_settings.cache_prefix / env_spec.kind
        pinfile_prefix = (
            self.workflow.deployment_settings.pinfile_prefix / env_spec.kind
        )
        deployment_prefix.mkdir(parents=True, exist_ok=True)
        cache_prefix.mkdir(parents=True, exist_ok=True)
        pinfile_prefix.mkdir(parents=True, exist_ok=True)

        env = env_spec.env_cls()(
            spec=env_spec,
            within=(
                self.get_env(
                    env_spec.within,
                    mountpoints=[deployment_prefix, cache_prefix, pinfile_prefix],
                )
                if env_spec.within is not None
                else None
            ),
            settings=self.workflow.software_deployment_provider_settings.get(
                self.plugins[env_spec.kind].name
            ),
            shell_executable=self.shell_executable,
            tempdir=Path(tempfile.gettempdir()),
            mountpoints=[self.workflow.source_cache_path, Path(os.getcwd())]
            + get_snakemake_searchpaths(),
            cache_prefix=cache_prefix,
            deployment_prefix=deployment_prefix,
            pinfile_prefix=pinfile_prefix,
        )
        if env in self.env_instances:
            # env with same content already instantiated, use that instead
            env = self.env_instances[env]

        self.specs_to_envs[env_spec] = env
        return env

    def register_in_global_variables(self, global_variables: Dict[str, Any]) -> None:
        for kind, plugin in self.plugins.items():
            if kind in global_variables:
                raise WorkflowError(
                    f"The software deployment plugin {plugin.name} provides an "
                    f"invalid kind {kind}, which is already reserved as a global "
                    "variable in snakemake."
                )

        for kind, plugin in self.plugins.items():

            def make_factory(plugin_val, kind_val):
                def factory(*args, within: EnvSpecBase, **kwargs):
                    if args:
                        raise WorkflowError(
                            "Positional arguments are not allowed in "
                            f"software deployment expressions (found in {kind}())"
                        )

                    for attr in plugin_val.env_spec_cls.source_path_attributes():
                        kwargs[attr] = EnvSpecSourceFile(kwargs[attr])

                    env_spec = plugin_val.env_spec_cls(*args, **kwargs)
                    env_spec.technical_init()

                    env_spec.within = within
                    env_spec.kind = kind_val
                    return env_spec

                return factory

            factory = make_factory(plugin, kind)

            args = list(get_function_params(plugin.env_spec_cls))
            args.append("within")
            overwrite_function_params(factory, args)
            global_variables[kind] = factory


@dataclass
class EnvSpecs:
    software_spec: Optional[Union[EnvSpecBase, Callable]] = None
    legacy_conda_env: Optional[Union[str, Path, Callable]] = None
    legacy_container_img: Optional[Union[str, Callable]] = None
    legacy_env_modules: Optional[Union[List[str], Callable]] = None

    def contains_conda_or_container(self) -> bool:
        def is_conda_or_container_spec(spec) -> bool:
            if spec is None:
                return False
            return (
                spec.kind in ("conda", "container")
                or is_conda_or_container_spec(spec.fallback)
                or is_conda_or_container_spec(spec.within)
            )

        return (
            is_conda_or_container_spec(self.software_spec)
            or self.legacy_conda_env is not None
            or self.legacy_container_img is not None
        )

    def is_empty(self) -> bool:
        return (
            self.software_spec is None
            and self.legacy_conda_env is None
            and self.legacy_container_img is None
            and self.legacy_env_modules is None
        )

    def is_callable(self) -> bool:
        return any(
            isinstance(spec, Callable)
            for spec in (
                self.software_spec,
                self.legacy_conda_env,
                self.legacy_container_img,
                self.legacy_env_modules,
            )
        )

    def resolve_callables(self, resolver: Callable) -> "EnvSpecs":
        resolved = copy(self)
        for spec in [
            "software_spec",
            "legacy_conda_env",
            "legacy_container_img",
            "legacy_env_modules",
        ]:
            if isinstance(getattr(self, spec), Callable):
                setattr(resolved, spec, resolver(getattr(self, spec)))

        return resolved

    def interpret(self) -> EnvSpecBase:
        assert not self.is_callable()

        if self.software_spec:
            # If a generic software spec is provided, always prefer that.
            return self.software_spec

        # Otherwise, use the old specs as fallback, converting to a generic spec.

        container_spec = (
            ContainerEnvSpec(self.legacy_container_img)
            if self.legacy_container_img is not None
            else None
        )
        conda_spec = None
        env_module_spec = None
        if self.legacy_conda_env is not None:
            spec = self.legacy_conda_env
            if isinstance(spec, Path):
                spec = str(spec)

            if spec.endswith(".yaml") or spec.endswith(".yml"):
                conda_spec = CondaEnvSpec(envfile=EnvSpecSourceFile(spec))
            elif is_local_file(spec) and os.path.isdir(spec):
                conda_spec = CondaEnvSpec(directory=Path(spec))
            else:
                conda_spec = CondaEnvSpec(name=spec)
            conda_spec.technical_init()

            # So far, if both container and conda were specified,
            # we assumed that the conda env is meant to be deployed
            # within the container, thus setting within below.
            conda_spec.within = container_spec

        if container_spec is not None:
            # If the container already contains the required software,
            # conda is just a fallback.
            container_spec.technical_init()
            container_spec.fallback = conda_spec

        if self.legacy_env_modules is not None:

            env_module_spec = EnvModuleEnvSpec(*self.legacy_env_modules)
            env_module_spec.technical_init()

            # If env modules shall be used, they are preferred (since they are infrastructure specific)
            # and the others are fallback.
            if conda_spec is not None:
                env_module_spec.fallback = conda_spec
            else:
                env_module_spec.fallback = container_spec

        if env_module_spec is not None:
            return env_module_spec
        if conda_spec is not None:
            return conda_spec
        if container_spec is not None:
            return container_spec
