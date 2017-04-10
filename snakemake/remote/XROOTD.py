__author__ = "Chris Burr"
__copyright__ = "Copyright 2017, Chris Burr"
__email__ = "christopher.burr@cern.ch"
__license__ = "MIT"

# built-ins
import os
from os.path import abspath, basename, dirname
import re
import http.client
import email.utils
# from itertools import product, chain
from contextlib import contextmanager

# module-specific
from snakemake.remote import AbstractRemoteProvider, DomainObject
from snakemake.exceptions import HTTPFileException, WorkflowError
import snakemake.io

import math
import email.utils
from time import mktime
import datetime
import functools
import concurrent.futures

# module-specific
from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import MissingOutputException, WorkflowError, WildcardError, RemoteFileException, XROOTDFileException

try:
    # third-party modules
    from XRootD import client
    from XRootD.client.flags import DirListFlags, OpenFlags, MkDirFlags, QueryCode
except ImportError as e:
    raise WorkflowError(
        "The Python 3 package 'XRootD' must be installed to use XROOTD "
        "remote() file functionality. %s" % e.msg
    )


class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, **kwargs):
        super(RemoteProvider, self).__init__(*args, **kwargs)

        self._xrd = XROOTDHelper(*args, **kwargs)

    def remote_interface(self):
        return self._xrd


class RemoteObject(AbstractRemoteObject):
    """ This is a class to interact with the AWS S3 object store.
    """

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)

        if provider:
            self._xrd = provider.remote_interface()
        else:
            self._xrd = XROOTDHelper(*args, **kwargs)

    # === Implementations of abstract class members ===

    def exists(self):
        return self._xrd.exists(self.file())

    def mtime(self):
        if self.exists():
            return self._xrd.file_last_modified(self.file())
        else:
            raise XROOTDFileException("The file does not seem to exist remotely: %s" % self.file())

    def size(self):
        if self.exists():
            return self._xrd.file_size(self.file())
        else:
            # TODO What is this????
            return self._iofile.size_local

    def download(self):
        self._xrd.download(self.file(), './'+self.file())

    def upload(self):
        self._xrd.upload('./'+self.file(), self.file())

    @property
    def list(self):
        return [f.name for f in self._xrd.list_directory(self.file())]


class XROOTDHelper(object):

    def __init__(self, url, port=1094):
        self._connection_string = 'root://{url}:{port}//'.format(url=url, port=port)
        self._client = client.FileSystem(self._connection_string)

    def exists(self, filename):
        status, _ = self._client.query(QueryCode.CHECKSUM, filename)
        if status.ok:
            return True
        elif status.errno == 3011:
            return False
        else:
            raise NotImplementedError('Got status object: '+repr(status))

    def _get_statinfo(self, filename):
        dirlist = self.list_directory(filename)
        matches = [f for f in dirlist.dirlist if f.name == basename(filename)]
        assert len(matches) == 1
        return matches[0].statinfo

    def file_last_modified(self, filename):
        # TODO Check this is in the right format
        return self._get_statinfo(filename).modtime

    def file_size(self, filename):
        # TODO Check this is in the right format (currently bytes)
        return self._get_statinfo(filename).size

    def _copy(self, source, destination):
        process = client.CopyProcess()
        process.add_job(self._connection_string+source, abspath(destination))
        process.prepare()
        status, returns = process.run()
        if not status.ok and returns[0]['status'].ok:
            raise XROOTDFileException('Error copying from '+source+' to '+destination)

    def download(self, source, destination):
        self._copy(self._connection_string+source, abspath(destination))

    def upload(self, source, destination):
        self._copy(abspath(destination), self._connection_string+source)

    def list_directory(self, filename):
        status, dirlist = self._client.dirlist(dirname(filename), DirListFlags.STAT)
        if not status.ok:
            raise XROOTDFileException('Error getting timestamp for '+filename, status, dirlist)
        return dirlist.dirlist
