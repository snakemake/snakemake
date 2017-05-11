__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import os
import sys
import re
from abc import ABCMeta, abstractmethod
from wrapt import ObjectProxy
import copy

# module-specific
import snakemake.io


class StaticRemoteObjectProxy(ObjectProxy):
    '''Proxy that implements static-ness for remote objects.

    The constructor takes a real RemoteObject and returns a proxy that
    behaves the same except for the exists() and mtime() methods.

    '''
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
    """ This is an abstract class to be used to derive remote provider classes. These might be used to hold common credentials,
        and are then passed to RemoteObjects.
    """
    __metaclass__ = ABCMeta

    def __init__(self, *args, keep_local=False, stay_on_remote=False, **kwargs):
        self.args = args
        self.stay_on_remote = stay_on_remote
        self.keep_local = keep_local
        self.kwargs = kwargs

    def remote(self, value, *args, keep_local=None, stay_on_remote=None, static=False, **kwargs):
        if snakemake.io.is_flagged(value, "temp"):
            raise SyntaxError(
                "Remote and temporary flags are mutually exclusive.")
        if snakemake.io.is_flagged(value, "protected"):
            raise SyntaxError(
                "Remote and protected flags are mutually exclusive.")
        if keep_local is None:
            keep_local = self.keep_local
        if stay_on_remote is None:
            stay_on_remote = self.stay_on_remote

        def _set_protocol(value):
            """Adds the default protocol to `value` if it doesn't already have one"""
            for protocol in self.available_protocols:
                if value.startswith(protocol):
                    break
            if value.startswith(protocol):
                value = value[len(protocol):]
                protocol = protocol
            else:
                protocol = self.default_protocol
            return protocol, value

        if isinstance(value, str):
            protocol, value = _set_protocol(value)
            value = protocol+value if stay_on_remote else value
        else:
            protocol, value = list(zip(*[_set_protocol(v) for v in value]))
            if len(set(protocol)) != 1:
                raise SyntaxError('A single protocol must be used per RemoteObject')
            protocol = set(protocol).pop()
            value = [protocol+v if stay_on_remote else v for v in value]

        provider = sys.modules[self.__module__]  # get module of derived class
        remote_object = provider.RemoteObject(
            *args, protocol=protocol, keep_local=keep_local, stay_on_remote=stay_on_remote,
            provider=provider.RemoteProvider(*self.args,  **self.kwargs), **kwargs
        )
        if static:
            remote_object = StaticRemoteObjectProxy(remote_object)
        return snakemake.io.flag(value, "remote_object", remote_object)

    def glob_wildcards(self, pattern, *args, **kwargs):
        args = self.args if not args else args
        kwargs = self.kwargs if not kwargs else kwargs

        referenceObj = snakemake.io.IOFile(self.remote(pattern, *args, **kwargs))

        if not referenceObj.remote_object.stay_on_remote:
            pattern = "./" + referenceObj.remote_object.name
            pattern = os.path.normpath(pattern)

        key_list = [k for k in referenceObj.remote_object.list]

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
    """ This is an abstract class to be used to derive remote object classes for
        different cloud storage providers. For example, there could be classes for interacting with
        Amazon AWS S3 and Google Cloud Storage, both derived from this common base class.
    """
    __metaclass__ = ABCMeta

    def __init__(self, *args, protocol=None, keep_local=False, stay_on_remote=False, provider=None, **kwargs):
        assert protocol is not None
        # self._iofile must be set before the remote object can be used, in io.py or elsewhere
        self._iofile = None
        self.args = args
        self.kwargs = kwargs

        self.keep_local = keep_local
        self.stay_on_remote = stay_on_remote
        self.provider = provider
        self.protocol = protocol

    @property
    def _file(self):
        if self._iofile is None:
            return None
        return self._iofile._file

    def file(self):
        return self._file

    def local_file(self):
        if self.stay_on_remote:
            return self._file[len(self.protocol):]
        else:
            return self._file

    def remote_file(self):
        return self.protocol+self.local_file()

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


class DomainObject(AbstractRemoteObject):
    """This is a mixin related to parsing components
        out of a location path specified as
        (host|IP):port/remote/location
    """
    def __init__(self, *args, **kwargs):
            super(DomainObject, self).__init__(*args, **kwargs)

    @property
    def _matched_address(self):
        return re.search("^(?P<host>[A-Za-z0-9\-\.]+)(?:\:(?P<port>[0-9]+))?(?P<path_remainder>.*)$", self.local_file())

    @property
    def name(self):
        return self.path_remainder

    @property
    def host(self):
        if self._matched_address:
            return self._matched_address.group("host")

    @property
    def port(self):
        return self._matched_address.group("port")

    @property
    def path_prefix(self):
        # this is the domain and port, however specified before the path remainder
        return self._iofile._file[:self._iofile._file.index(self.path_remainder)]

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
