__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os
import re
import ftplib
import collections
from itertools import chain
from contextlib import contextmanager

# module-specific
from snakemake.remote import AbstractRemoteProvider, DomainObject
from snakemake.exceptions import FTPFileException, WorkflowError

try:
    # third-party modules
    import ftputil
    import ftputil.session
except ImportError as e:
    raise WorkflowError("The Python 3 package 'ftputil' " +
        "must be installed to use SFTP remote() file functionality. %s" % e.msg)


class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, stay_on_remote=False, immediate_close=False, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)

        self.immediate_close = immediate_close

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return 'ftp://'

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ['ftp://', 'ftps://']

    def remote(self, value, *args, encrypt_data_channel=None, immediate_close=None, **kwargs):
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, collections.Iterable):
            values = value
        else:
            raise TypeError('Invalid type ({}) passed to remote: {}'.format(type(value), value))

        for i, file in enumerate(values):
            match = re.match('^(ftps?)://.+', file)
            if match:
                protocol, = match.groups()
                if protocol == 'ftps' and encrypt_data_channel:
                    raise SyntaxError('encrypt_data_channel=False cannot be used with a ftps:// url')
                if protocol == 'ftp' and encrypt_data_channel not in [None, False]:
                    raise SyntaxError('encrypt_data_channel=Trie cannot be used with a ftp:// url')
            else:
                if encrypt_data_channel:
                    values[i] = 'ftps://' + file
                else:
                    values[i] = 'ftp://' + file

        should_close = immediate_close if immediate_close else self.immediate_close
        return super(RemoteProvider, self).remote(values, *args, encrypt_data_channel=encrypt_data_channel, immediate_close=should_close, **kwargs)


class RemoteObject(DomainObject):
    """ This is a class to interact with an FTP server.
    """

    def __init__(self, *args, keep_local=False, provider=None, encrypt_data_channel=False, immediate_close=False, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)

        self.encrypt_data_channel = encrypt_data_channel
        self.immediate_close      = immediate_close

    def close(self):
        if hasattr(self, "conn") and isinstance(self.conn, ftputil.FTPHost) and not self.immediate_close:
            try:
                self.conn.keep_alive()
                self.conn.close()
            except:
                pass

    # === Implementations of abstract class members ===

    @contextmanager #makes this a context manager. after 'yield' is __exit__()
    def ftpc(self):
        if (not hasattr(self, "conn") or (hasattr(self, "conn") and not isinstance(self.conn, ftputil.FTPHost))) or self.immediate_close:
            # if args have been provided to remote(), use them over those given to RemoteProvider()
            args_to_use = self.provider.args
            if len(self.args):
                args_to_use = self.args

            # use kwargs passed in to remote() to override those given to the RemoteProvider()
            # default to the host and port given as part of the file, falling back to one specified
            # as a kwarg to remote() or the RemoteProvider (overriding the latter with the former if both)
            kwargs_to_use = {}
            kwargs_to_use["host"] = self.host
            kwargs_to_use["username"] = None
            kwargs_to_use["password"] = None
            kwargs_to_use["port"] = int(self.port) if self.port else 21
            kwargs_to_use["encrypt_data_channel"] = self.encrypt_data_channel

            for k,v in self.provider.kwargs.items():
                kwargs_to_use[k] = v
            for k,v in self.kwargs.items():
                kwargs_to_use[k] = v

            ftp_base_class = ftplib.FTP_TLS if kwargs_to_use["encrypt_data_channel"] else ftplib.FTP

            ftp_session_factory = ftputil.session.session_factory(
                           base_class=ftp_base_class,
                           port=kwargs_to_use["port"],
                           encrypt_data_channel= kwargs_to_use["encrypt_data_channel"],
                           debug_level=None)

            conn = ftputil.FTPHost(kwargs_to_use["host"], kwargs_to_use["username"], kwargs_to_use["password"], session_factory=ftp_session_factory)
            if self.immediate_close:
                yield conn
            else:    
                self.conn = conn
                yield self.conn
        elif not self.immediate_close:
            yield self.conn
        
        # after returning from the context manager, close the connection if the scope is local
        if self.immediate_close:
            try:
                conn.keep_alive()
                conn.close()
            except:
                pass

    def exists(self):
        if self._matched_address:
            with self.ftpc() as ftpc:
                return ftpc.path.exists(self.remote_path)
            return False
        else:
            raise FTPFileException("The file cannot be parsed as an FTP path in form 'host:port/abs/path/to/file': %s" % self.local_file())

    def mtime(self):
        if self.exists():
            with self.ftpc() as ftpc:
                try:
                    # requires write access
                    ftpc.synchronize_times()
                except:
                    pass
                return ftpc.path.getmtime(self.remote_path)
        else:
            raise FTPFileException("The file does not seem to exist remotely: %s" % self.local_file())

    def size(self):
        if self.exists():
            with self.ftpc() as ftpc:
                return ftpc.path.getsize(self.remote_path)
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        with self.ftpc() as ftpc:
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
                os.sync() # ensure flush to disk
            else:
                raise FTPFileException("The file does not seem to exist remotely: %s" % self.local_file())

    def upload(self):
        with self.ftpc() as ftpc:
            ftpc.synchronize_times()
            ftpc.upload(source=self.local_path, target=self.remote_path)

    @property
    def list(self):
        file_list = []

        first_wildcard = self._iofile.constant_prefix()
        dirname = first_wildcard.replace(self.path_prefix, "")

        with self.ftpc() as ftpc:
            file_list = [(os.path.join(dirpath, f) if dirpath != "." else f)
                    for dirpath, dirnames, filenames in ftpc.walk(dirname)
                    for f in chain(filenames, dirnames)]
            file_list = [file_path[1:] if file_path[0] == "/" else file_path for file_path in file_list]

        return file_list
