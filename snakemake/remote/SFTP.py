__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os

# module-specific
from snakemake.remote import AbstractRemoteProvider, PooledDomainObject
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


class RemoteObject(PooledDomainObject):
    """This is a class to interact with an SFTP server."""

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )

    def get_default_kwargs(self, **defaults):
        """define defaults beyond those set in PooledDomainObject"""
        return super().get_default_kwargs(
            **{
                "port": 22,
            }
        )

    def create_connection(self, *args_to_use, **kwargs_to_use):
        """open an SFTP connection"""
        return pysftp.Connection(*args_to_use, **kwargs_to_use)

    def close_connection(self, connection):
        """close an SFTP connection"""
        connection.close()

    # === Implementations of abstract class members ===

    def exists(self):
        if self._matched_address:
            with self.connection_pool.item() as sftpc:
                return sftpc.exists(self.remote_path)
            return False
        else:
            raise SFTPFileException(
                "The file cannot be parsed as an SFTP path in form 'host:port/path/to/file': %s"
                % self.local_file()
            )

    def mtime(self):
        if self.exists():
            with self.connection_pool.item() as sftpc:
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
        with self.connection_pool.item() as sftpc:
            return (
                sftpc.stat(self.remote_path).st_mtime > time
                or sftpc.lstat(self.remote_path).st_mtime > time
            )

    def size(self):
        if self.exists():
            with self.connection_pool.item() as sftpc:
                attr = sftpc.stat(self.remote_path)
                return int(attr.st_size)
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        with self.connection_pool.item() as sftpc:
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

        with self.connection_pool.item() as sftpc:
            for part in list_remote_dir:
                try:
                    sftpc.chdir(part)
                except IOError:
                    sftpc.mkdir(part)
                    sftpc.chdir(part)

    def upload(self):
        if self.provider.mkdir_remote:
            self.mkdir_remote_path()

        with self.connection_pool.item() as sftpc:
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

        with self.connection_pool.item() as sftpc:

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
