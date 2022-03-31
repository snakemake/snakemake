__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2022, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import collections
from itertools import chain
from os import path

import snakemake
from snakemake.common import parse_uri
from snakemake.common.future import cached_property
from snakemake.common.plugin import (
    PluginException,
    find_plugins,
    internal_submodules,
    load_plugins,
    verify_plugin,
)
from snakemake.common.remote import (
    AbstractRemoteObject,
    AbstractRemoteProvider,
    AbstractRemoteRetryObject,
    DomainObject,
    PooledDomainObject,
    StaticRemoteObjectProxy,
    check_deprecated_retry,
)
from snakemake.logging import logger
from snakemake.utils import min_version, raise_


def snakemake_min_version(mod):
    """raises WorkFlowError"""
    min_version(mod.__min_snakemake_version__)


plugin_checks = {
    snakemake_min_version: None,
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

remote_plugin_prefix = "snakemake-plugin-remote-"


def verify_remote_plugin(mod):
    verify_plugin(
        mod,
        extra_attrs=["RemoteProvider", "RemoteObject"],
        checks=plugin_checks,
        prefix=remote_plugin_prefix,
    )


# Load remote plugins in a list,
# Also have them available for direct importing
# by setting them on the globals() dict.
plugin_remote_modules = list(
    load_plugins(
        plugin_modules=find_plugins(prefix="snakemake-plugin-remote-"),
        globals_dict=globals(),
        verify_func=verify_remote_plugin,
    )
)


class AutoRemoteProvider:
    """Automatically gather all RemoteProviders and find the correct one based on the protocol"""

    @cached_property
    def protocol_mapping(self):
        remote_path_list = [path.join(p, "remote") for p in snakemake.__path__]
        provider_list = [
            module.RemoteProvider
            for module in chain(
                internal_submodules(remote_path_list), plugin_remote_modules
            )
        ]

        # assemble scheme mapping
        mapping = {}
        for Provider in provider_list:
            try:
                available_protocols = Provider.available_protocols.__get__(Provider)
                for protocol in available_protocols:
                    if protocol[-3:] != "://":
                        logger.warning(
                            f"RemoteProvider from {Provider.__module__} "
                            "has a protocol {protocol} that does not end with ://."
                        )
                        continue
                    protocol_short = protocol[:-3]  # remove "://" suffix
                    mapping[protocol_short] = Provider
            except Exception:
                logger.debug(f"Something went wrong with loading {Provider}")
                continue
        return maping

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
