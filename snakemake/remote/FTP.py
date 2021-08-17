__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os
import re
import ftplib
import collections
from itertools import chain

# module-specific
from snakemake.remote import AbstractRemoteProvider, PooledDomainObject

from snakemake.exceptions import FTPFileException, WorkflowError
from snakemake.utils import os_sync

try:
    # third-party modules
    import ftputil
    import ftputil.session
except ImportError as e:
    raise WorkflowError(
        "The Python 3 package 'ftputil' "
        + "must be installed to use FTP remote() file functionality. %s" % e.msg
    )


class RemoteProvider(AbstractRemoteProvider):

    supports_default = True
    allows_directories = True

    def __init__(
        self, *args, keep_local=False, stay_on_remote=False, is_default=False, **kwargs
    ):
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs
        )

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "ftp://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["ftp://", "ftps://"]

    def remote(self, value, *args, encrypt_data_channel=None, **kwargs):
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, collections.abc.Iterable):
            values = value
        else:
            raise TypeError(
                "Invalid type ({}) passed to remote: {}".format(type(value), value)
            )

        for i, file in enumerate(values):
            match = re.match("^(ftps?)://.+", file)
            if match:
                (protocol,) = match.groups()
                if protocol == "ftps" and encrypt_data_channel in [None, False]:
                    raise SyntaxError(
                        "encrypt_data_channel=False cannot be used with a ftps:// url"
                    )
                if protocol == "ftp" and encrypt_data_channel not in [None, False]:
                    raise SyntaxError(
                        "encrypt_data_channel=True cannot be used with a ftp:// url"
                    )
            else:
                if encrypt_data_channel:
                    values[i] = "ftps://" + file
                else:
                    values[i] = "ftp://" + file

        values = [
            super(RemoteProvider, self).remote(
                value, *args, encrypt_data_channel=encrypt_data_channel, **kwargs
            )
            for value in values
        ]
        if len(values) == 1:
            return values[0]
        else:
            return values


class RemoteObject(PooledDomainObject):
    """This is a class to interact with an FTP server."""

    def __init__(
        self,
        *args,
        keep_local=False,
        provider=None,
        encrypt_data_channel=False,
        **kwargs
    ):
        # just save to kwargs, but leave in function def so it's explicit
        kwargs["encrypt_data_channel"] = encrypt_data_channel
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )

    # === Implementations of abstract class members ===

    def get_default_kwargs(self, **defaults):
        """define defaults beyond thos set in PooledDomainObject"""
        return super().get_default_kwargs(
            **{"port": 21, "password": None, "username": None}
        )

    @property
    def conn_keywords(self):
        """returns list of keywords relevant to a unique connection"""
        return ["host", "port", "username", "encrypt_data_channel"]

    def create_connection(self, *args_to_use, **kwargs_to_use):
        """create a connection to the FTP server using the saved arguments"""
        ftp_base_class = (
            ftplib.FTP_TLS if kwargs_to_use["encrypt_data_channel"] else ftplib.FTP
        )

        ftp_session_factory = ftputil.session.session_factory(
            base_class=ftp_base_class,
            port=kwargs_to_use["port"],
            encrypt_data_channel=kwargs_to_use["encrypt_data_channel"],
            debug_level=None,
        )

        return ftputil.FTPHost(
            kwargs_to_use["host"],
            kwargs_to_use["username"],
            kwargs_to_use["password"],
            session_factory=ftp_session_factory,
        )

    def close_connection(self, connection):
        """close an FTP connection"""
        connection.keep_alive()
        connection.close()

    def exists(self):
        if self._matched_address:
            with self.connection_pool.item() as ftpc:
                return ftpc.path.exists(self.remote_path)
            return False
        else:
            raise FTPFileException(
                "The file cannot be parsed as an FTP path in form 'host:port/abs/path/to/file': %s"
                % self.local_file()
            )

    def mtime(self):
        if self.exists():
            with self.connection_pool.item() as ftpc:
                try:
                    # requires write access
                    ftpc.synchronize_times()
                except:
                    pass
                return ftpc.path.getmtime(self.remote_path)
        else:
            raise FTPFileException(
                "The file does not seem to exist remotely: %s" % self.local_file()
            )

    def size(self):
        if self.exists():
            with self.connection_pool.item() as ftpc:
                return ftpc.path.getsize(self.remote_path)
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        with self.connection_pool.item() as ftpc:
            if self.exists():
                # if the destination path does not exist
                if make_dest_dirs:
                    os.makedirs(os.path.dirname(self.local_path), exist_ok=True)
                try:
                    # requires write access
                    ftpc.synchronize_times()
                except:
                    pass
                ftpc.download(source=self.remote_path, target=self.local_path)
                os_sync()  # ensure flush to disk
            else:
                raise FTPFileException(
                    "The file does not seem to exist remotely: %s" % self.local_file()
                )

    def upload(self):
        with self.connection_pool.item() as ftpc:
            ftpc.synchronize_times()
            ftpc.upload(source=self.local_path, target=self.remote_path)

    @property
    def list(self):
        file_list = []

        first_wildcard = self._iofile.constant_prefix()
        dirname = first_wildcard.replace(self.path_prefix, "")

        with self.connection_pool.item() as ftpc:
            file_list = [
                (os.path.join(dirpath, f) if dirpath != "." else f)
                for dirpath, dirnames, filenames in ftpc.walk(dirname)
                for f in chain(filenames, dirnames)
            ]
            file_list = [
                file_path[1:] if file_path[0] == "/" else file_path
                for file_path in file_list
            ]

        return file_list
