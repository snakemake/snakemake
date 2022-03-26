__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2022, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

from os import path
import collections
from snakemake.common import parse_uri
from snakemake.logging import logger
from snakemake.common.plugin import find_plugins, verify_module_attrs, PluginException, internal_submodules
from .common import AbstractRemoteObject, AbstractRemoteProvider, AbstractRemoteRetryObject, check_deprecated_retry, DomainObject, PooledDomainObject, StaticRemoteObjectProxy
import snakemake
from itertools import chain


plugin_remote_modules = []

for plugin_module in find_plugins(prefix='snakemake-plugin-remote-'):
    try:
        verify_module_attrs(
            plugin_module, extra_attrs=["RemoteProvider", "RemoteOjbect"]
        )
    except PluginException as e:
        logger.warning(f"Plugin {plugin_module.name} incorrect: {e}")
        continue

    if not isinstance(plugin_module.RemoteObject, AbstractRemoteObject):
        raise PluginException(
            f"RemoteObject from plugin {plugin_module.name} "
            "is not an instance of AbstractRemoteObject"
        )

    if not isinstance(plugin_module.RemoteProvider, AbstractRemoteProvider):
        raise PluginException(
            f"RemoteProvider from plugin {plugin_module.name} "
            "is not an instance of AbstractRemoteProvider"
        )

    if plugin_module.name in globals():
        logger.debug("Plugin {plugin_module.name} already loaded.")
        continue

    plugin_remote_modules.append(plugin_module)
    globals()[plugin_module.name] = plugin_module
    globals()[plugin_module.name.lower()] = plugin_module


class AutoRemoteProvider:
    def __init__(self):
        """Automatically gather all RemoteProviders"""
        remote_path_list = [path.join(p, "remote") for p in snakemake.__path__]
        provider_list = [
            module.RemoteProvider for module in
            chain(internal_submodules(remote_path_list), plugin_remote_modules)
        ]

        # assemble scheme mapping
        self.protocol_mapping = {}
        for Provider in provider_list:
            for protocol in Provider().available_protocols:
                if protocol[-3:] != '://':
                    logger.warning(
                        f"RemoteProvider from {Provider.__module__} "
                        "has a protocol {protocol} that does not end with ://."
                    )
                    continue
                protocol_short = protocol[:-3]  # remove "://" suffix
                self.protocol_mapping[protocol_short] = Provider

    def remote(self, value, *args, provider_kws=None, **kwargs):
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, collections.abc.Iterable):
            values = value
        else:
            raise TypeError(
                f"Invalid type ({type(value)}) passed to remote: {value}"
            )

        provider_remote_list = []
        for value in values:
            # select provider
            o = parse_uri(value)
            Provider = self.protocol_mapping.get(o.scheme)

            if Provider is None:
                raise TypeError(f"Could not find remote provider for: {value}")

            # use provider's remote
            provider_kws = {} if provider_kws is None else provider_kws.copy()

            provider_remote_list.append(
                Provider(**provider_kws).remote(value, *args, **kwargs)
            )

        return (
            provider_remote_list[0]
            if len(provider_remote_list) == 1
            else provider_remote_list
        )


AUTO = AutoRemoteProvider()
