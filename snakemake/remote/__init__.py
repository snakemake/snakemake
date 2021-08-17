__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import os
import sys
import re
from functools import partial
from abc import ABCMeta, abstractmethod
from wrapt import ObjectProxy
from contextlib import contextmanager

try:
    from connection_pool import ConnectionPool
except ImportError:
    # we just won't pool connections if it's not installed
    #  Should there be a warning? Should there be a runtime flag?
    pass
import copy
import collections

# module-specific
import snakemake.io
from snakemake.logging import logger
from snakemake.common import parse_uri


class StaticRemoteObjectProxy(ObjectProxy):
    """Proxy that implements static-ness for remote objects.

    The constructor takes a real RemoteObject and returns a proxy that
    behaves the same except for the exists() and mtime() methods.

    """

    def exists(self):
        return True

    def mtime(self):
        return float("-inf")

    def is_newer(self, time):
        return False

    def __copy__(self):
        copied_wrapped = copy.copy(self.__wrapped__)
        return type(self)(copied_wrapped)

    def __deepcopy__(self):
        copied_wrapped = copy.deepcopy(self.__wrapped__)
        return type(self)(copied_wrapped)


class AbstractRemoteProvider:
    """This is an abstract class to be used to derive remote provider classes. These might be used to hold common credentials,
    and are then passed to RemoteObjects.
    """

    __metaclass__ = ABCMeta

    supports_default = False
    allows_directories = False

    def __init__(
        self, *args, keep_local=False, stay_on_remote=False, is_default=False, **kwargs
    ):
        self.args = args
        self.stay_on_remote = stay_on_remote
        self.keep_local = keep_local
        self.is_default = is_default
        self.kwargs = kwargs

    def remote(
        self, value, *args, keep_local=None, stay_on_remote=None, static=False, **kwargs
    ):
        if snakemake.io.is_flagged(value, "temp"):
            raise SyntaxError("Remote and temporary flags are mutually exclusive.")
        if snakemake.io.is_flagged(value, "protected"):
            raise SyntaxError("Remote and protected flags are mutually exclusive.")
        if keep_local is None:
            keep_local = self.keep_local
        if stay_on_remote is None:
            stay_on_remote = self.stay_on_remote

        def _set_protocol(value):
            """Adds the default protocol to `value` if it doesn't already have one"""
            protocol = self.default_protocol
            for p in self.available_protocols:
                if value.startswith(p):
                    value = value[len(p) :]
                    protocol = p
                    break
            return protocol, value

        if isinstance(value, str):
            protocol, value = _set_protocol(value)
            value = protocol + value if stay_on_remote else value
        else:
            protocol, value = list(zip(*[_set_protocol(v) for v in value]))
            if len(set(protocol)) != 1:
                raise SyntaxError("A single protocol must be used per RemoteObject")
            protocol = set(protocol).pop()
            value = [protocol + v if stay_on_remote else v for v in value]

        if "protocol" not in kwargs:
            if "protocol" not in self.kwargs:
                kwargs["protocol"] = protocol
            else:
                kwargs["protocol"] = self.kwargs["protocol"]

        provider = sys.modules[self.__module__]  # get module of derived class
        remote_object = provider.RemoteObject(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            provider=self,
            **kwargs,
        )
        if static:
            remote_object = StaticRemoteObjectProxy(remote_object)
        return snakemake.io.flag(value, "remote_object", remote_object)

    def glob_wildcards(self, pattern, *args, **kwargs):
        args = self.args if not args else args
        kwargs = self.kwargs if not kwargs else kwargs

        referenceObj = snakemake.io._IOFile(self.remote(pattern, *args, **kwargs))
        remote_object = snakemake.io.get_flag_value(referenceObj, "remote_object")
        if not remote_object.stay_on_remote:
            pattern = "./" + remote_object.name
            pattern = os.path.normpath(pattern)

        key_list = [k for k in remote_object.list]

        return snakemake.io.glob_wildcards(pattern, files=key_list)

    @abstractmethod
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        pass

    @abstractmethod
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        pass

    @abstractmethod
    def remote_interface(self):
        pass


class AbstractRemoteObject:
    """This is an abstract class to be used to derive remote object classes for
    different cloud storage providers. For example, there could be classes for interacting with
    Amazon AWS S3 and Google Cloud Storage, both derived from this common base class.
    """

    __metaclass__ = ABCMeta

    def __init__(
        self,
        *args,
        protocol=None,
        keep_local=False,
        stay_on_remote=False,
        provider=None,
        **kwargs,
    ):
        assert protocol is not None
        # self._iofile must be set before the remote object can be used, in io.py or elsewhere
        self._iofile = None
        self.args = args
        self.kwargs = kwargs

        self.keep_local = keep_local
        self.stay_on_remote = stay_on_remote
        self.provider = provider
        self.protocol = protocol

    async def inventory(self, cache: snakemake.io.IOCache):
        """From this file, try to find as much existence and modification date
        information as possible.
        """
        # If this is implemented in a remote object, results have to be stored in
        # the given IOCache object.
        pass

    @abstractmethod
    def get_inventory_parent(self):
        pass

    @property
    def _file(self):
        if self._iofile is None:
            return None
        return self._iofile._file

    def file(self):
        return self._file

    def local_file(self):
        if self.stay_on_remote:
            return self._file[len(self.protocol) :]
        else:
            return self._file

    def remote_file(self):
        return self.protocol + self.local_file()

    @abstractmethod
    def close(self):
        pass

    @abstractmethod
    def exists(self):
        pass

    @abstractmethod
    def mtime(self):
        pass

    @abstractmethod
    def size(self):
        pass

    @abstractmethod
    def download(self, *args, **kwargs):
        pass

    @abstractmethod
    def upload(self, *args, **kwargs):
        pass

    @abstractmethod
    def list(self, *args, **kwargs):
        pass

    @abstractmethod
    def name(self, *args, **kwargs):
        pass

    @abstractmethod
    def remote(self, value, keep_local=False, stay_on_remote=False):
        pass

    @abstractmethod
    def remove(self):
        raise NotImplementedError("Removal of files is unavailable for this remote")

    def local_touch_or_create(self):
        self._iofile.touch_or_create()


class DomainObject(AbstractRemoteObject):
    """This is a mixin related to parsing components
    out of a location path specified as
    (host|IP):port/remote/location
    """

    def __init__(self, *args, **kwargs):
        super(DomainObject, self).__init__(*args, **kwargs)

    @property
    def _matched_address(self):
        return re.search(
            r"^(?P<protocol>[a-zA-Z]+\://)?(?P<host>[A-Za-z0-9\-\.]+)(?:\:(?P<port>[0-9]+))?(?P<path_remainder>.*)$",
            self.local_file(),
        )

    @property
    def name(self):
        return self.path_remainder

    # if we ever parse out the protocol directly
    # @property
    # def protocol(self):
    #    if self._matched_address:
    #        return self._matched_address.group("protocol")

    @property
    def host(self):
        if self._matched_address:
            return self._matched_address.group("host")

    @property
    def port(self):
        if self._matched_address:
            return self._matched_address.group("port")

    @property
    def path_prefix(self):
        # this is the domain and port, however specified before the path remainder
        return self._iofile._file[: self._iofile._file.index(self.path_remainder)]

    @property
    def path_remainder(self):
        if self._matched_address:
            return self._matched_address.group("path_remainder")

    @property
    def local_path(self):
        return self._iofile._file

    @property
    def remote_path(self):
        return self.path_remainder


class PooledDomainObject(DomainObject):
    """This adds conection pooling to DomainObjects
    out of a location path specified as
    (host|IP):port/remote/location
    """

    connection_pools = {}

    def __init__(self, *args, pool_size=100, immediate_close=False, **kwargs):
        super(PooledDomainObject, self).__init__(*args, **kwargs)
        self.pool_size = 100
        self.immediate_close = immediate_close

    def get_default_kwargs(self, **defaults):
        defaults.setdefault("host", self.host)
        defaults.setdefault("port", int(self.port) if self.port else None)
        return defaults

    def get_args_to_use(self):
        """merge the objects args with the parent provider

        Positional Args: use those of object or fall back to ones from provider
        Keyword Args: merge with any defaults
        """
        # if args have been provided to remote(),
        #  use them over those given to RemoteProvider()
        args_to_use = self.provider.args
        if len(self.args):
            args_to_use = self.args

        # use kwargs passed in to remote() to override those given to the RemoteProvider()
        #  default to the host and port given as part of the file,
        #  falling back to one specified as a kwarg to remote() or the RemoteProvider
        #  (overriding the latter with the former if both)
        kwargs_to_use = self.get_default_kwargs()
        for k, v in self.provider.kwargs.items():
            kwargs_to_use[k] = v
        for k, v in self.kwargs.items():
            kwargs_to_use[k] = v

        return args_to_use, kwargs_to_use

    @contextmanager
    def get_connection(self):
        """get a connection from a pool or create a new one"""
        if not self.immediate_close and "connection_pool" in sys.modules:
            # if we can (and the user doesn't override) use a pool
            with self.connection_pool.item() as conn:
                yield conn
        else:
            # otherwise create a one-time connection
            args_to_use, kwargs_to_use = self.get_args_to_use()
            conn = self.create_connection(*args_to_use, **kwargs_to_use)
            try:
                yield conn
            finally:
                self.close_connection(conn)

    @property
    def conn_keywords(self):
        """returns list of keywords relevant to a unique connection"""
        return ["host", "port", "username"]

    @property
    def connection_pool(self):
        """set up a pool of re-usable active connections"""
        # merge this object's values with those of its parent provider
        args_to_use, kwargs_to_use = self.get_args_to_use()

        # hashing connection pool on tuple of relevant arguments. There
        # may be a better way to do this
        conn_pool_label_tuple = (
            type(self),
            *args_to_use,
            *[kwargs_to_use.get(k, None) for k in self.conn_keywords],
        )

        if conn_pool_label_tuple not in self.connection_pools:
            create_callback = partial(
                self.create_connection, *args_to_use, **kwargs_to_use
            )
            self.connection_pools[conn_pool_label_tuple] = ConnectionPool(
                create_callback, close=self.close_connection, max_size=self.pool_size
            )

        return self.connection_pools[conn_pool_label_tuple]

    @abstractmethod
    def create_connection(self):
        """handle the protocol specific job of creating a connection"""
        pass

    @abstractmethod
    def close_connection(self, connection):
        """handle the protocol specific job of closing a connection"""
        pass


class AutoRemoteProvider:
    @property
    def protocol_mapping(self):
        # automatically gather all RemoteProviders
        import pkgutil
        import importlib.util

        provider_list = []
        for remote_submodule in pkgutil.iter_modules(snakemake.remote.__path__):
            path = (
                os.path.join(remote_submodule.module_finder.path, remote_submodule.name)
                + ".py"
            )

            module_name = remote_submodule.name
            spec = importlib.util.spec_from_file_location(module_name, path)
            module = importlib.util.module_from_spec(spec)

            try:
                sys.modules[module_name] = module
                spec.loader.exec_module(module)
            except Exception as e:
                logger.debug(f"Autoloading {module_name} failed: {e}")
                continue

            provider_list.append(module.RemoteProvider)

        # assemble scheme mapping
        protocol_dict = {}
        for Provider in provider_list:
            for protocol in Provider().available_protocols:
                protocol_short = protocol[:-3]  # remove "://" suffix
                protocol_dict[protocol_short] = Provider

        return protocol_dict

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
