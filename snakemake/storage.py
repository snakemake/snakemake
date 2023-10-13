import copy, sys
import subprocess
from dataclasses import dataclass, field
from typing import Any, Dict, Optional
from snakemake.workflow import Workflow
from snakemake_interface_common.exceptions import WorkflowError, InvalidPluginException
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry
from snakemake_interface_storage_plugins.storage_provider import StorageProviderBase
from snakemake_interface_storage_plugins.storage_object import (
    StorageObjectWrite,
    StorageObjectRead,
)
from snakemake_interface_executor_plugins.settings import DeploymentMethod
from snakemake.common import __version__


class StorageRegistry:
    attrs = {
        "workflow",
        "_storages",
        "_default_storage_provider",
        "default_storage_provider",
        "_register_default_storage",
        "register_storage",
        "infer_provider",
        "_storage_object",
        "__getattribute__",
        "attrs",
    }

    def __init__(self, workflow: Workflow):
        self.workflow = workflow
        self._storages = dict()
        self._default_storage_provider = None

        if self.workflow.storage_settings.default_storage_provider is not None:
            self._register_default_storage()

    def _register_default_storage(self):
        plugin_name = self.workflow.storage_settings.default_storage_provider
        if (
            not StoragePluginRegistry().is_installed(plugin_name)
            and self.workflow.deployment_settings.default_storage_provider_auto_deploy
        ):
            try:
                subprocess.run(
                    ["pip", "install", f"snakemake-storage-plugin-{plugin_name}"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            except subprocess.CalledProcessError as e:
                raise WorkflowError(
                    f"Failed to install storage plugin {plugin_name} via pip: {e.stdout.decode()}",
                    e,
                )
            StoragePluginRegistry().collect_plugins()
        self._default_storage_provider = self.register_storage(
            plugin_name,
            is_default=True,
        )

    @property
    def default_storage_provider(self):
        return self._default_storage_provider

    def register_storage(
        self,
        provider: Optional[str] = None,
        tag: Optional[str] = None,
        is_default: bool = False,
        **settings,
    ):
        if provider is None:
            raise WorkflowError('Storage provider must be specified (provider="...").')
        if provider != provider.lower():
            raise WorkflowError("Storage provider must be lowercase.")

        # First retrieve plugin in order to ensure that it is there.
        plugin = StoragePluginRegistry().get_plugin(provider)

        tagged_settings = self.workflow.storage_provider_settings.get(provider)
        if tagged_settings is None:
            final_settings = None
        else:
            final_settings = tagged_settings.get_settings(tag)
        if final_settings is None:
            final_settings = plugin.settings_cls()

        final_settings = copy.copy(final_settings)
        final_settings.__dict__.update(**settings)

        name = tag if tag else plugin.name

        local_prefix = self.workflow.storage_settings.local_storage_prefix / name

        if is_default and not (
            issubclass(
                plugin.storage_provider.get_storage_object_cls(), StorageObjectWrite
            )
            and issubclass(
                plugin.storage_provider.get_storage_object_cls(), StorageObjectRead
            )
        ):
            raise WorkflowError(
                "Default storage provider must be a read-write storage provider, but "
                f"{plugin.name} is not."
            )

        provider_instance = plugin.storage_provider(
            local_prefix=local_prefix,
            settings=final_settings,
            keep_local=self.workflow.storage_settings.keep_storage_local,
            is_default=is_default,
        )
        self._storages[name] = provider_instance
        return provider_instance

    def infer_provider(self, query: str):
        plugins = [
            plugin.name
            for plugin in StoragePluginRegistry().plugins.values()
            if plugin.storage_provider.is_valid_query(query)
        ]
        if len(plugins) == 1:
            return plugins[0]
        elif len(plugins) == 0:
            raise WorkflowError(
                f"No storage provider found for query {query}. "
                "Either install the required storage plugin or check your query. "
                "Also consider to explictly specify the storage provider to get a more "
                "informative error message."
            )
        else:
            raise WorkflowError(
                f"Multiple suitable storage providers found for query {query}: {', '.join(plugins)}. "
                "Explictly specify the storage provider."
            )

    def __getattribute__(self, name: str) -> Any:
        if name == "attrs":
            return super().__getattribute__(name)
        elif name in self.attrs:
            return super().__getattribute__(name)
        else:
            return StorageProviderProxy(registry=self, name=name)

    def __call__(
        self,
        query: str,
        retrieve: bool = True,
        keep_local: bool = False,
    ):
        return self._storage_object(
            query, provider=None, retrieve=retrieve, keep_local=keep_local
        )

    def _storage_object(
        self,
        query: str,
        provider: Optional[str] = None,
        retrieve: bool = True,
        keep_local: bool = False,
    ):
        provider_name = provider

        if provider_name is None:
            provider_name = self.infer_provider(query)

        provider = self._storages.get(provider_name)
        if provider is None:
            provider = self.register_storage(provider_name)

        return provider.object(query, retrieve=retrieve, keep_local=keep_local)


class StorageProviderProxy:
    def __init__(self, registry: StorageRegistry, name: str):
        self.registry = registry
        self.name = name

    def __call__(
        self,
        query: str,
        retrieve: bool = True,
        keep_local: bool = False,
    ):
        return self.registry._storage_object(
            query, provider=self.name, retrieve=retrieve, keep_local=keep_local
        )
