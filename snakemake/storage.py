import copy
from dataclasses import dataclass, field
from typing import Any, Dict, Optional
from snakemake.workflow import Workflow
from snakemake_interface_common.exceptions import WorkflowError
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry
from snakemake_interface_storage_plugins.storage_provider import StorageProviderBase
from snakemake_interface_storage_plugins.storage_object import StorageObjectReadWrite


@dataclass
class StorageRegistry:
    workflow: Workflow
    _storages: Dict[str, StorageProviderBase] = field(default_factory=dict)
    _default_storage_provider: Optional[StorageProviderBase] = None

    def __post_init__(self):
        if self.workflow.storage_settings.default_storage_provider is not None:
            self._default_storage_provider = self.register_storage(
                self.workflow.storage_settings.default_storage_provider, is_default=True
            )

    @property
    def default_storage_provider(self):
        return self._default_storage_provider

    def register_storage(
        self,
        provider: str,
        tag: Optional[str] = None,
        is_default: bool = False,
        **settings,
    ):
        if provider != provider.lower():
            raise WorkflowError("Storage provider must be lowercase.")

        # First retrieve plugin in order to ensure that it is there.
        plugin = StoragePluginRegistry().get_plugin(provider)

        final_settings = self.workflow.storage_provider_settings.get(provider)
        if final_settings is None:
            final_settings = plugin.settings_cls()
        else:
            final_settings = final_settings.get_settings(tag)

        final_settings = copy.copy(final_settings)
        final_settings.__dict__.update(**settings)

        name = f"{plugin.name}:{tag}"

        local_prefix = self.workflow.storage_settings.local_storage_prefix / name

        if is_default and not isinstance(
            plugin.storage_provider.storage_object_cls, StorageObjectReadWrite
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

    def storage_object(
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
            if ":" in provider_name:
                raise WorkflowError(
                    "Usage of tagged storage (e.g. s3:mytag) requires a storage "
                    "directive (storage: ...) to be used before in the workflow that "
                    "defines the tag and the associated storage provider and settings."
                )
            provider = self.register_storage(provider_name)

        return provider.object(query, retrieve=retrieve, keep_local=keep_local)
