__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import os, sys
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

    def remote(self, value, keep_local=False, additional_args=None, additional_kwargs=None):
        additional_args = [] if not additional_args else additional_args
        additional_kwargs = {} if not additional_kwargs else additional_kwargs

        if snakemake.io.is_flagged(value, "temp"):
            raise SyntaxError(
                "Remote and temporary flags are mutually exclusive.")
        if snakemake.io.is_flagged(value, "protected"):
            raise SyntaxError(
                "Remote and protected flags are mutually exclusive.")

        provider = sys.modules[self.__module__] # get module of derived class
        remote_object = provider.RemoteObject(*additional_args, keep_local=keep_local, provider=provider.RemoteProvider(*self.args,  **self.kwargs), **additional_kwargs)

        return snakemake.io.flag(
                value, 
                "remote_object",
                remote_object
            )

    def glob_wildcards(self, pattern, additional_args=None, additional_kwargs=None):
        additional_args   = self.args if not additional_args else additional_args
        additional_kwargs = self.kwargs if not additional_kwargs else additional_kwargs
        
        referenceObj = snakemake.io.IOFile(self.remote(pattern, additional_args, additional_kwargs))

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
