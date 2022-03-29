__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2022, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

from os import path
import collections
from snakemake.common import parse_uri
from snakemake.logging import logger
from snakemake.common.plugin import (
    find_plugins,
    PluginException,
    internal_submodules,
    load_plugins,
)
from .common import (
    AbstractRemoteObject,
    AbstractRemoteProvider,
    AbstractRemoteRetryObject,
    check_deprecated_retry,
    DomainObject,
    PooledDomainObject,
    StaticRemoteObjectProxy,
)
import snakemake
from itertools import chain
from snakemake.utils import min_version, raise_

plugin_checks = {
    lambda mod: min_version(mod.__min_snakemake_version__): None,
    lambda mod: not issubclass(
        mod.RemoteObject, AbstractRemoteObject
    ): lambda mod: raise_(
        PluginException(
            "RemoteObject from plugin {mod.__name__} "
            "is not an instance of AbstractRemoteObject"
        )
    ),
    lambda mod: not issubclass(
        mod.RemoteProvider, AbstractRemoteProvider
    ): lambda mod: raise_(
        PluginException(
            f"RemoteProvider from plugin {mod.__name__} "
            "is not an instance of AbstractRemoteProvider"
        )
    ),
}
# Load remote plugins in a list,
# Also have them available for direct importing
# by setting them on the globals() dict.
plugin_remote_modules = list(
    load_plugins(
        plugin_modules=find_plugins(prefix="snakemake-plugin-remote-"),
        globals_dict=globals(),
        extra_attrs=["RemoteProvider", "RemoteOjbect", "__min_snakemake_version__"],
        checks=plugin_checks,
    )
)


class AutoRemoteProvider:
    def __init__(self):
        """Automatically gather all RemoteProviders"""
        remote_path_list = [path.join(p, "remote") for p in snakemake.__path__]
        provider_list = [
            module.RemoteProvider
            for module in chain(
                internal_submodules(remote_path_list), plugin_remote_modules
            )
        ]

        # assemble scheme mapping
        self.protocol_mapping = {}
        for Provider in provider_list:
            for protocol in Provider().available_protocols:
                if protocol[-3:] != "://":
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
            raise TypeError(f"Invalid type ({type(value)}) passed to remote: {value}")

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
