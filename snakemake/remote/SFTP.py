__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os, re
from contextlib import contextmanager

# module-specific
from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import SFTPFileException
import snakemake.io 

try:
    # third-party modules
    import pysftp
except ImportError as e:
    raise WorkflowError("The Python 3 package 'pysftp' " + 
        "must be installed to use SFTP remote() file functionality. %s" % e.msg)


class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, **kwargs):
        super(RemoteProvider, self).__init__(*args, **kwargs)

class RemoteObject(AbstractRemoteObject):
    """ This is a class to interact with the AWS S3 object store.
    """

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)
        
    # === Implementations of abstract class members ===

    @contextmanager #makes this a context manager. after 'yield' is __exit__()
    def sftpc(self):     
        # if args have been provided to remote(), use them over those given to RemoteProvider()
        args_to_use = self.provider.args
        if len(self.args):
            args_to_use = self.args

        # use kwargs passed in to remote() to override those given to the RemoteProvider()
        # default to the host and port given as part of the file, falling back to one specified
        # as a kwarg to remote() or the RemoteProvider (overriding the latter with the former if both)
        kwargs_to_use = {}
        kwargs_to_use["host"] = self.host
        kwargs_to_use["port"] = self.port
        for k,v in self.provider.kwargs.items():
            kwargs_to_use[k] = v
        for k,v in self.kwargs.items():
            kwargs_to_use[k] = v

        conn = pysftp.Connection(*args_to_use, **kwargs_to_use)
        yield conn
        conn.close()

    def exists(self):
        if self._matched_address:
            with self.sftpc() as sftpc:
                return sftpc.exists(self.remote_path)
                if sftpc.exists(self.remote_path):
                    return sftpc.isfile(self.remote_path)
            return False
        else:
            raise SFTPFileException("The file cannot be parsed as an s3 path in form 'host:port/path/to/file': %s" % self.file())

    def mtime(self):
        if self.exists():
            with self.sftpc() as sftpc:
                attr = sftpc.stat(self.remote_path)
                return int(attr.st_mtime)
        else:
            raise SFTPFileException("The file does not seem to exist remotely: %s" % self.file())

    def size(self):
        if self.exists():
            with self.sftpc() as sftpc:
                attr = sftpc.stat(self.remote_path)
                return int(attr.st_size)
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        with self.sftpc() as sftpc:
            if self.exists():
                # if the destination path does not exist
                if not os.path.exists(os.path.dirname(self.local_path)) and make_dest_dirs:
                    os.makedirs(os.path.dirname(self.local_path))

                sftpc.get(remotepath=self.remote_path, localpath=self.local_path, preserve_mtime=True)
            else:
                raise SFTPFileException("The file does not seem to exist remotely: %s" % self.file())
    def upload(self):
        with self.sftpc() as sftpc:
            sftpc.put(localpath=self.local_path, remotepath=self.remote_path, confirm=True, preserve_mtime=True)

    @property
    def list(self):
        file_list = []

        pattern = os.path.normpath(self.remote_path)
        first_wildcard = re.search("{[^{]", pattern)
        dirname = os.path.dirname(pattern[:first_wildcard.start(
        )]) if first_wildcard else os.path.dirname(pattern)
        if not dirname:
            dirname = "."

        with self.sftpc() as sftpc:
            def _append_item(file_path):
                file_path = file_path[1:] if file_path[0] == "/" else file_path
                file_list.append(file_path)

            sftpc.walktree(dirname, fcallback=_append_item, dcallback=_append_item, ucallback=_append_item)

        return file_list

    # === Related methods ===

    @property
    def name(self):
        return self.path_remainder

    @property
    def _matched_address(self):
        return re.search("^(?:(?P<protocol>.*)\://)?(?P<host>[A-Za-z0-9\-\.]+)(?:\:(?P<port>[0-9]+))?(?P<path_remainder>.*)$", self._iofile._file)
    
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
        if self._matched_address.group("port"):
            return int(self._matched_address.group("port"))
        else:
            return 22
    
    @property
    def path_remainder(self):
        if self._matched_address:
            return self._matched_address.group("path_remainder")

    @property
    def local_path(self):
        return self._iofile._file
        #return "{host}{path_remainder}".format(host=self.host, port=self.port, path_remainder=self.path_remainder)    

    @property
    def remote_path(self):
        return self.path_remainder

