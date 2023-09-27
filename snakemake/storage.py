import copy
from dataclasses import dataclass, field
from typing import Any, Dict, Optional
from snakemake.workflow import Workflow
from snakemake_interface_common.plugin_registry.plugin import TaggedSettings
from snakemake_interface_common.exceptions import WorkflowError
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry
from snakemake_interface_storage_plugins.storage_provider import StorageProviderBase


@dataclass
class StorageRegistry:
    workflow: Workflow
    _storages: Dict[str, StorageProviderBase] = field(default_factory=dict)
    _default_storage_provider: Optional[StorageProviderBase] = None

    def __post_init__(self):
        if self.workflow.storage_settings.default_storage_provider is not None:
            self._default_storage_provider = self.register_storage(self.workflow.storage_settings.default_storage_provider)
    
    @property
    def default_storage_provider(self):
        return self._default_storage_provider

    def register_storage(self, provider: str, tag: Optional[str] = None, **settings):
        if not provider.is_lower():
            raise WorkflowError("Storage provider must be lowercase.")
        
        # First retrieve plugin in order to ensure that it is there.
        plugin = StoragePluginRegistry().get_plugin(provider)
        final_settings = copy.copy(self.workflow.storage_provider_settings[provider].get_settings(tag))
        final_settings.__dict__.update(**settings)

        name = f"{plugin.name}:{tag}"

        provider_instance = plugin.storage_provider(settings=settings)
        self._storages[name] = provider_instance
        return provider_instance

    def storage_object(
        self,
        query: Any,
        storage_name: str,
        retrieve: bool = True,
        keep_local: bool = False,
    ):
        provider = self._storages.get(storage_name)

        if provider is None:
            raise WorkflowError(
                f"Storage {storage_name} is not registered. Please use the 'storage' "
                "directive to register this storage."
            )
        
        return provider.object(query, retrieve=retrieve, keep_local=keep_local)