__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os, sys
from abc import ABCMeta, abstractmethod

import snakemake.io

class RemoteProvider:
    """ This is an abstract class to be used to derive remote provider classes. These might be used to hold common credentials,
        and are then passed to RemoteObjects.
    """
    __metaclass__ = ABCMeta

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def remote(self, value, provider=None, keep_local=False, additional_args=None, additional_kwargs=None):
        # set some defaults
        additional_args   = self.args if not additional_args else additional_args
        additional_kwargs = self.kwargs if not additional_kwargs else additional_kwargs
        provider          = sys.modules[self.__module__] if not provider else provider

        return snakemake.io.remote(value=value, provider=provider, keep_local=keep_local, additional_args=additional_args, additional_kwargs=additional_kwargs)

    def glob_wildcards(self, pattern, provider=None, additional_args=None, additional_kwargs=None):
        additional_args   = self.args if not additional_args else additional_args
        additional_kwargs = self.kwargs if not additional_kwargs else additional_kwargs
        
        referenceObj = snakemake.io.IOFile(self.remote(pattern, additional_args, additional_kwargs))

        pattern = "./"+ referenceObj._remote_object.name
        pattern = os.path.normpath(pattern)

        key_list = [k.name for k in referenceObj._remote_object.list] 

        return snakemake.io.glob_wildcards(pattern, files=key_list)

class RemoteObject:
    """ This is an abstract class to be used to derive remote object classes for 
        different cloud storage providers. For example, there could be classes for interacting with 
        Amazon AWS S3 and Google Cloud Storage, both derived from this common base class.
    """
    __metaclass__ = ABCMeta

    def __init__(self, ioFile, *args, **kwargs):
        self._iofile = ioFile
        self._file = ioFile._file

    @abstractmethod
    def file(self):
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
    def remote(self, value, keep_local=False):
        pass