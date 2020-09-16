__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os
from contextlib import contextmanager

# module-specific
from snakemake.remote import AbstractRemoteProvider, DomainObject
from snakemake.exceptions import SFTPFileException, WorkflowError
from snakemake.utils import os_sync

try:
    # third-party modules
    import pysftp
except ImportError as e:
    raise WorkflowError(
        "The Python 3 package 'pysftp' "
        + "must be installed to use SFTP remote() file functionality. %s" % e.msg
    )


class RemoteProvider(AbstractRemoteProvider):

    supports_default = True
    allows_directories = True

    def __init__(
        self,
        *args,
        keep_local=False,
        stay_on_remote=False,
        is_default=False,
        mkdir_remote=True,
        **kwargs
    ):
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs
        )

        self.mkdir_remote = mkdir_remote

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "sftp://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["ssh://", "sftp://"]


class RemoteObject(DomainObject):
    """This is a class to interact with an SFTP server."""

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )

    # === Implementations of abstract class members ===

    @contextmanager  # makes this a context manager. after 'yield' is __exit__()
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
        kwargs_to_use["port"] = int(self.port) if self.port else 22
        for k, v in self.provider.kwargs.items():
            kwargs_to_use[k] = v
        for k, v in self.kwargs.items():
            kwargs_to_use[k] = v

        conn = pysftp.Connection(*args_to_use, **kwargs_to_use)
        yield conn
        conn.close()

    def exists(self):
        if self._matched_address:
            with self.sftpc() as sftpc:
                return sftpc.exists(self.remote_path)
            return False
        else:
            raise SFTPFileException(
                "The file cannot be parsed as an SFTP path in form 'host:port/path/to/file': %s"
                % self.local_file()
            )

    def mtime(self):
        if self.exists():
            with self.sftpc() as sftpc:
                # As per local operation, don't follow symlinks when reporting mtime
                attr = sftpc.lstat(self.remote_path)
                return int(attr.st_mtime)
        else:
            raise SFTPFileException(
                "The file does not seem to exist remotely: %s" % self.local_file()
            )

    def is_newer(self, time):
        """Returns true if the file is newer than time, or if it is
        a symlink that points to a file newer than time."""
        with self.sftpc() as sftpc:
            return (
                sftpc.stat(self.remote_path).st_mtime > time
                or sftpc.lstat(self.remote_path).st_mtime > time
            )

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
                if make_dest_dirs:
                    os.makedirs(os.path.dirname(self.local_path), exist_ok=True)

                sftpc.get(
                    remotepath=self.remote_path,
                    localpath=self.local_path,
                    preserve_mtime=True,
                )
                os_sync()  # ensure flush to disk
            else:
                raise SFTPFileException(
                    "The file does not seem to exist remotely: %s" % self.local_file()
                )

    def mkdir_remote_path(self):
        remote_dir = os.path.dirname(self.remote_path)
        list_remote_dir = []
        while True:
            remote_dir, base = os.path.split(remote_dir)
            if not base and remote_dir:
                list_remote_dir.insert(0, remote_dir)
                break
            list_remote_dir.insert(0, base)

        with self.sftpc() as sftpc:
            for part in list_remote_dir:
                try:
                    sftpc.chdir(part)
                except IOError:
                    sftpc.mkdir(part)
                    sftpc.chdir(part)

    def upload(self):
        if self.provider.mkdir_remote:
            self.mkdir_remote_path()

        with self.sftpc() as sftpc:
            sftpc.put(
                localpath=self.local_path,
                remotepath=self.remote_path,
                confirm=True,
                preserve_mtime=True,
            )

    @property
    def list(self):
        file_list = []

        first_wildcard = self._iofile.constant_prefix()
        dirname = first_wildcard.replace(self.path_prefix, "")

        with self.sftpc() as sftpc:

            def _append_item(file_path):
                file_path = file_path.lstrip("/")
                file_list.append(file_path)

            sftpc.walktree(
                dirname,
                fcallback=_append_item,
                dcallback=_append_item,
                ucallback=_append_item,
            )

        return file_list
