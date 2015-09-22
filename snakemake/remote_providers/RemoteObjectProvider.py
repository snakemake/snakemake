__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

from abc import ABCMeta, abstractmethod


class RemoteObject:
    """ This is an abstract class to be used to derive remote object classes for 
        different cloud storage providers. For example, there could be classes for interacting with 
        Amazon AWS S3 and Google Cloud Storage, both derived from this common base class.
    """
    __metaclass__ = ABCMeta

    def __init__(self, ioFile):
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
