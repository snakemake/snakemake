__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import re

from snakemake.remote_providers.RemoteObjectProvider import RemoteObject
from snakemake.exceptions import MissingOutputException, WorkflowError, WildcardError, RemoteFileException, S3FileException
from snakemake.remote_providers.implementations.S3 import S3Helper

import boto


class RemoteObject(RemoteObject):
    """ This is a class to interact with the AWS S3 object store.
    """

    def __init__(self, *args, **kwargs):
        super(RemoteObject, self).__init__(*args, **kwargs)

        # pass all args but the first, which is the ioFile
        self._s3c = S3Helper(*args[1:], **kwargs)

    # === Implementations of abstract class members ===

    def file(self):
        return self._file

    def exists(self):
        if self._matched_s3_path:
            return self._s3c.exists_in_bucket(self.s3_bucket, self.s3_key)
        else:
            raise S3FileException("The file cannot be parsed as an s3 path in form 'bucket/key': %s" % self.file())

    def mtime(self):
        if self.exists():
            return self._s3c.key_last_modified(self.s3_bucket, self.s3_key)
        else:
            raise S3FileException("The file does not seem to exist remotely: %s" % self.file())

    def size(self):
        if self.exists():
            return self._s3c.key_size(self.s3_bucket, self.s3_key)
        else:
            return self._iofile.size_local

    def download(self):
        self._s3c.download_from_s3(self.s3_bucket, self.s3_key, self.file())

    def upload(self):
        conn = boto.connect_s3()
        if self.size() > 5000:
            self._s3c.upload_to_s3_multipart(self.s3_bucket, self.file(), self.s3_key)
        else:
            self._s3c.upload_to_s3(self.s3_bucket, self.file(), self.s3_key)

    @property
    def list(self):
        return self._s3c.list_keys(self.s3_bucket)

    # === Related methods ===

    @property
    def _matched_s3_path(self):
        return re.search("(?P<bucket>[^/]*)/(?P<key>.*)", self.file())

    @property
    def s3_bucket(self):
        if len(self._matched_s3_path.groups()) == 2:
            return self._matched_s3_path.group("bucket")
        return None

    @property
    def name(self):
        return self.s3_key

    @property
    def s3_key(self):
        if len(self._matched_s3_path.groups()) == 2:
            return self._matched_s3_path.group("key")

    def s3_create_stub(self):
        if self._matched_s3_path:
            if not self.exists:
                self._s3c.download_from_s3(self.s3_bucket, self.s3_key, self.file, createStubOnly=True)
        else:
            raise S3FileException("The file to be downloaded cannot be parsed as an s3 path in form 'bucket/key': %s" %
                                  self.file())
