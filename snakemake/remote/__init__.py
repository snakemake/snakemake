__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import os, sys, re
from abc import ABCMeta, abstractmethod

# module-specific
import snakemake.io
from snakemake.exceptions import RemoteFileException

class AbstractRemoteProvider:
    """ This is an abstract class to be used to derive remote provider classes. These might be used to hold common credentials,
        and are then passed to RemoteObjects.
    """
    __metaclass__ = ABCMeta

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def remote(self, value, *args, keep_local=False, **kwargs):
        if snakemake.io.is_flagged(value, "temp"):
            raise SyntaxError(
                "Remote and temporary flags are mutually exclusive.")
        if snakemake.io.is_flagged(value, "protected"):
            raise SyntaxError(
                "Remote and protected flags are mutually exclusive.")

        provider = sys.modules[self.__module__] # get module of derived class
        remote_object = provider.RemoteObject(*args, keep_local=keep_local, provider=provider.RemoteProvider(*self.args,  **self.kwargs), **kwargs)

        return snakemake.io.flag(
                value, 
                "remote_object",
                remote_object
            )

    def glob_wildcards(self, pattern, *args, **kwargs):
        args   = self.args if not args else args
        kwargs = self.kwargs if not kwargs else kwargs
        
        referenceObj = snakemake.io.IOFile(self.remote(pattern, *args, **kwargs))

        pattern = "./"+ referenceObj.remote_object.name
        pattern = os.path.normpath(pattern)

        key_list = [k for k in referenceObj.remote_object.list] 

        return snakemake.io.glob_wildcards(pattern, files=key_list)

    @abstractmethod
    def remote_interface(self):
        pass

class AbstractRemoteObject:
    """ This is an abstract class to be used to derive remote object classes for 
        different cloud storage providers. For example, there could be classes for interacting with 
        Amazon AWS S3 and Google Cloud Storage, both derived from this common base class.
    """
    __metaclass__ = ABCMeta

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        # self._iofile must be set before the remote object can be used, in io.py or elsewhere
        self._iofile = None
        self.args = args
        self.kwargs = kwargs

        self.keep_local = keep_local
        self.provider = provider

    @property
    def _file(self):
        return self._iofile._file
    
    def file(self):
        return self._file

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
    def remote(self, value, keep_local=False):
        pass

class DomainObject(AbstractRemoteObject):
    """This is a mixin related to parsing components
        out of a location path specified as
        (host|IP):port/remote/location
    """
    def __init__(self, *args, **kwargs):
            super(DomainObject, self).__init__(*args, **kwargs)

    @property
    def _matched_address(self):
        return re.search("^(?P<host>[A-Za-z0-9\-\.]+)(?:\:(?P<port>[0-9]+))?(?P<path_remainder>.*)$", self._iofile._file)

    @property
    def name(self):
        return self.path_remainder
    
    @property
    def protocol(self):
        if self._matched_address:
            return self._matched_address.group("protocol")

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
