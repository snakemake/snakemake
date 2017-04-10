__author__ = "Chris Burr"
__copyright__ = "Copyright 2017, Chris Burr"
__email__ = "christopher.burr@cern.ch"
__license__ = "MIT"

import os
from os.path import abspath, basename, dirname, normpath, isdir

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError, XROOTDFileException

try:
    # third-party modules
    from XRootD import client
    from XRootD.client.flags import DirListFlags, MkDirFlags
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
    """ This is a class to interact with XROOTD servers.
    """

    def __init__(self, *args, keep_local=False, use_remote=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, use_remote=use_remote, provider=provider, **kwargs)

        if provider:
            self._xrd = provider.remote_interface()
        else:
            self._xrd = XROOTDHelper(*args, **kwargs)

    @property
    def _file(self):
        if self._iofile is None:
            return None
        return normpath('./'+self._iofile._file)

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
        return basename(filename) in [f.name for f in self.list_directory(dirname(filename)+'/')]

    def _get_statinfo(self, filename):
        matches = [f for f in self.list_directory('/'+filename) if f.name == basename(filename)]
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
        process.add_job(source, destination)
        process.prepare()
        status, returns = process.run()
        if not status.ok or not returns[0]['status'].ok:
            raise XROOTDFileException('Error copying from '+source+' to '+destination, repr(status), repr(returns))

    def download(self, source, destination):
        dest_dir = abspath(dirname(destination))
        if not isdir(dest_dir):
            os.makedirs(dest_dir)
        self._copy(self._connection_string+source, abspath(destination))

    def makedirs(self, dirname):
        status, _ = self._client.mkdir('/'+dirname+'/', MkDirFlags.MAKEPATH)
        if not status.ok:
            print('*'*100, repr(status))
            raise XROOTDFileException('Failed to create directory /'+dirname, repr(status))

    def upload(self, source, destination):
        dest_dir = normpath(dirname(destination))
        # TODO Make this more like isdir
        if not self.exists(dest_dir):
            self.makedirs(dest_dir)
        self._copy(abspath(source), self._connection_string+destination)

    def list_directory(self, filename):
        status, dirlist = self._client.dirlist(dirname('/'+filename), DirListFlags.STAT)
        if not status.ok:
            raise XROOTDFileException('Error getting timestamp for '+filename, status, dirlist)
        return dirlist.dirlist
