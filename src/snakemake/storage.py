import copy
from typing import Any, List, Mapping, Optional, Union
from snakemake.io import flag
from snakemake.workflow import Workflow
from snakemake_interface_common.exceptions import WorkflowError
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry
from snakemake_interface_storage_plugins.storage_provider import StorageProviderBase
from snakemake_interface_storage_plugins.storage_object import (
    StorageObjectWrite,
    StorageObjectRead,
)
from snakemake.io import MaybeAnnotated
from snakemake.common import __version__
from snakemake.logging import logger


def flag_with_storage_object(path: MaybeAnnotated, storage_object):
    modified = flag(storage_object.local_path(), "storage_object", storage_object)
    modified.flags.update(getattr(path, "flags", {}).copy())

    return modified


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
        self._storages: Mapping[str, StorageProviderBase] = dict()
        self._default_storage_provider = None

        if self.workflow.storage_settings.default_storage_provider is not None:
            self._default_storage_provider = self.register_storage(
                self.workflow.storage_settings.default_storage_provider,
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
        if final_settings is None and plugin.settings_cls is not None:
            final_settings = plugin.settings_cls()

        if final_settings is not None:
            final_settings = copy.copy(final_settings)
            final_settings.__dict__.update(**settings)

            plugin.validate_settings(final_settings)

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

        keep_local = settings.get(
            "keep_local", self.workflow.storage_settings.keep_storage_local
        )
        retrieve = settings.get(
            "retrieve", self.workflow.storage_settings.retrieve_storage
        )
        provider_instance = plugin.storage_provider(
            logger=logger,
            local_prefix=local_prefix,
            settings=final_settings,
            keep_local=keep_local,
            retrieve=retrieve,
            is_default=is_default,
            wait_for_free_local_storage=self.workflow.storage_settings.wait_for_free_local_storage,
        )
        self._storages[name] = provider_instance
        # if a tagged storage provider is registered before the untagged then the
        # untagged provider is registered later without the settings
        # prevent the settings loss by registering it here
        if tag is not None and plugin.name not in self._storages:
            local_prefix = (
                self.workflow.storage_settings.local_storage_prefix / plugin.name
            )
            provider_instance = plugin.storage_provider(
                logger=logger,
                local_prefix=local_prefix,
                settings=final_settings,
                keep_local=keep_local,
                retrieve=retrieve,
                is_default=is_default,
                wait_for_free_local_storage=self.workflow.storage_settings.wait_for_free_local_storage,
            )
            self._storages[plugin.name] = provider_instance
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
                "Also consider to explicitly specify the storage provider to get a more "
                "informative error message."
            )
        else:
            raise WorkflowError(
                f"Multiple suitable storage providers found for query {query}: {', '.join(plugins)}. "
                "Explicitly specify the storage provider."
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
        retrieve: Optional[bool] = None,
        keep_local: Optional[bool] = None,
    ):
        return self._storage_object(
            query, provider=None, retrieve=retrieve, keep_local=keep_local
        )

    def _storage_object(
        self,
        query: Union[str, List[str]],
        provider: Optional[str] = None,
        retrieve: Optional[bool] = None,
        keep_local: Optional[bool] = None,
    ):
        if isinstance(query, list):
            return [
                self._storage_object(
                    q, provider=provider, retrieve=retrieve, keep_local=keep_local
                )
                for q in query
            ]

        provider_name = provider

        if provider_name is None:
            provider_name = self.infer_provider(query)

        provider = self._storages.get(provider_name)
        if provider is None:
            provider = self.register_storage(provider_name)

        query_validity = provider.is_valid_query(query)
        if not query_validity:
            raise WorkflowError(
                f"Error applying storage provider {provider_name} "
                "(see https://snakemake.github.io/snakemake-plugin-catalog/plugins/"
                f"storage/{provider}.html). {query_validity}"
            )

        storage_object = provider.object(
            query, retrieve=retrieve, keep_local=keep_local
        )

        return flag_with_storage_object(storage_object.query, storage_object)


class StorageProviderProxy:
    def __init__(self, registry: StorageRegistry, name: str):
        self.registry = registry
        self.name = name

    def __call__(
        self,
        query: str,
        retrieve: Optional[bool] = None,
        keep_local: Optional[bool] = None,
    ):
        return self.registry._storage_object(
            query, provider=self.name, retrieve=retrieve, keep_local=keep_local
        )
