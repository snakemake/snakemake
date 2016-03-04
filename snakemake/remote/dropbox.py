__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os, re
from contextlib import contextmanager

# module-specific
from snakemake.remote import AbstractRemoteProvider, AbstractRemoteObject
from snakemake.exceptions import DropboxFileException, WorkflowError
import snakemake.io

try:
    # third-party modules
    import dropbox # The official Dropbox API library
except ImportError as e:
    raise WorkflowError("The Python 3 package 'dropbox' "
                        "must be installed to use Dropbox remote() file "
                        "functionality. %s" % e.msg)


class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, **kwargs):
        super(RemoteProvider, self).__init__(*args, **kwargs)

        self._dropboxc = dropbox.Dropbox(*args, **kwargs)
        try:
            self._dropboxc.users_get_current_account()
        except dropbox.exceptions.AuthError as err:
                DropboxFileException("ERROR: Invalid Dropbox OAuth access token; try re-generating an access token from the app console on the web.")

    def remote_interface(self):
        return self._dropboxc

class RemoteObject(AbstractRemoteObject):
    """ This is a class to interact with the AWS S3 object store.
    """

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)

        if provider:
            self._dropboxc = provider.remote_interface()
        else:
            self._dropboxc = dropbox.Dropbox(*args, **kwargs)
            try:
                self._dropboxc.users_get_current_account()
            except dropbox.exceptions.AuthError as err:
                    DropboxFileException("ERROR: Invalid Dropbox OAuth access token; try re-generating an access token from the app console on the web.")

    # === Implementations of abstract class members ===

    def exists(self):
        try:
            metadata = self._dropboxc.files_get_metadata(self.remote_file())
            return True
        except:
            return False

    def mtime(self):
        if self.exists():
            metadata = self._dropboxc.files_get_metadata(self.remote_file())
            epochTime = metadata.server_modified.timestamp()
            return epochTime
        else:
            raise DropboxFileException("The file does not seem to exist remotely: %s" % self.remote_file())

    def size(self):
        if self.exists():
            metadata = self._dropboxc.files_get_metadata(self.remote_file())
            return int(metadata.size)
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        if self.exists():
            # if the destination path does not exist, make it
            if make_dest_dirs:
                os.makedirs(os.path.dirname(self.file()), exist_ok=True)

            self._dropboxc.files_download_to_file(self.file(), self.remote_file())
        else:
            raise DropboxFileException("The file does not seem to exist remotely: %s" % self.remote_file())

    def upload(self, mode=dropbox.files.WriteMode('overwrite')):
        size = os.path.getsize(self.file())
        # Chunk file into 10MB slices because Dropbox does not accept more than 150MB chunks
        chunksize = 10000000
        with open(self.file(), mode='rb') as f:
            data = f.read(chunksize)
            # Start upload session
            res = self._dropboxc.files_upload_session_start(data)
            offset = len(data)

            # Upload further chunks until file is complete
            while len(data) == chunksize:
                data = f.read(chunksize)
                self._dropboxc.files_upload_session_append(data, res.session_id, offset)
                offset += len(data)

            # Finish session and store in the desired path
            self._dropboxc.files_upload_session_finish(
                f.read(chunksize),
                dropbox.files.UploadSessionCursor(res.session_id, offset),
                dropbox.files.CommitInfo(path=self.remote_file(), mode=mode))

    def remote_file(self):
        return "/"+self.file() if not self.file().startswith("/") else self.file()

    @property
    def name(self):
        return self.file()

    @property
    def list(self):
        file_list = []

        first_wildcard = self._iofile.constant_prefix()
        dirname = "/" + first_wildcard if not first_wildcard.startswith("/") else first_wildcard

        while '//' in dirname:
            dirname = dirname.replace('//', '/')
        dirname = dirname.rstrip('/')

        for item in self._dropboxc.files_list_folder(dirname, recursive=True).entries:
            file_list.append( os.path.join(os.path.dirname(item.path_lower), item.name).lstrip("/") )

        return file_list
