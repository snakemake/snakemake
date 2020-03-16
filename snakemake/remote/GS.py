__author__ = "Johannes Köster"
__copyright__ = "Copyright 2017-2019, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import base64
import os
import re
import struct

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError
from snakemake.common import lazy_property

try:
    import google.cloud
    from google.cloud import storage
    from crc32c import crc32
except ImportError as e:
    raise WorkflowError(
        "The Python 3 packages 'google-cloud-sdk' and `crc32c` "
        "need to be installed to use GS remote() file functionality. %s" % e.msg
    )


class Crc32cCalculator:
    """The Google Python client doesn't provide a way to stream a file being
       written, so we can wrap the file object in an additional class to
       do custom handling. This is so we don't need to download the file
       and then stream read it again to calculate the hash.
   """

    def __init__(self, fileobj):
        self._fileobj = fileobj
        self.digest = 0

    def write(self, chunk):
        self._fileobj.write(chunk)
        self._update(chunk)

    def _update(self, chunk):
        """Given a chunk from the read in file, update the hexdigest
        """
        self.digest = crc32(chunk, self.digest)

    def hexdigest(self):
        """Return the hexdigest of the hasher.
           The Base64 encoded CRC32c is in big-endian byte order.
           See https://cloud.google.com/storage/docs/hashes-etags
        """
        return base64.b64encode(struct.pack(">I", self.digest)).decode("utf-8")


class RemoteProvider(AbstractRemoteProvider):

    supports_default = True

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

        self.client = storage.Client(*args, **kwargs)

    def remote_interface(self):
        return self.client

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "gs://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["gs://"]


class RemoteObject(AbstractRemoteObject):
    def __init__(
        self, *args, keep_local=False, provider=None, user_project=None, **kwargs
    ):
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )

        if provider:
            self.client = provider.remote_interface()
        else:
            self.client = storage.Client(*args, **kwargs)

        # keep user_project available for when bucket is initialized
        self._user_project = user_project

        self._key = None
        self._bucket_name = None
        self._bucket = None
        self._blob = None

    # === Implementations of abstract class members ===

    def exists(self):
        return self.blob.exists()

    def mtime(self):
        if self.exists():
            self.update_blob()
            t = self.blob.updated
            return t.timestamp()
        else:
            raise WorkflowError(
                "The file does not seem to exist remotely: %s" % self.local_file()
            )

    def size(self):
        if self.exists():
            self.update_blob()
            return self.blob.size // 1024
        else:
            return self._iofile.size_local

    def download(self, retry_count=3):
        if not self.exists():
            return None

        # Create the directory for the intended file
        os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)
        checksums_match = False

        # Use boolean to not stress filesystem on each loop
        while not checksums_match:

            # ideally we could calculate hash while streaming to file with provided function
            # https://github.com/googleapis/python-storage/issues/29
            with open(self.local_file(), "wb") as blob_file:
                parser = Crc32cCalculator(blob_file)
                self.blob.download_to_file(parser)
            os.sync()

            # Compute local hash and verify correct
            if parser.hexdigest() != self.blob.crc32c:
                os.remove(self.local_file())
                retry_count -= 1
            else:
                checksums_match = True

            # After three tries, likely not going to work
            if retry_count == 0:
                raise WorkflowError(
                    "Failed to download {} from remote: checksums did not match {} times.".format(
                        self.remote_file(), retry_count
                    )
                )

        return self.local_file()

    def upload(self):
        try:
            if not self.bucket.exists():
                self.bucket.create()
                self.update_blob()
            self.blob.upload_from_filename(self.local_file())
        except google.cloud.exceptions.Forbidden as e:
            raise WorkflowError(
                e,
                "When running locally, make sure that you are authenticated "
                "via gcloud (see Snakemake documentation). When running in a "
                "kubernetes cluster, make sure that storage-rw is added to "
                "--scopes (see Snakemake documentation).",
            )

    @property
    def name(self):
        return self.key

    @property
    def list(self):
        return [k.name for k in self.bucket.list_blobs()]

    # ========= Helpers ===============

    def update_blob(self):
        self._blob = self.bucket.get_blob(self.key)

    @lazy_property
    def bucket(self):
        return self.client.bucket(self.bucket_name, user_project=self._user_project)

    @lazy_property
    def blob(self):
        return self.bucket.blob(self.key)

    @lazy_property
    def bucket_name(self):
        return self.parse().group("bucket")

    @lazy_property
    def key(self):
        return self.parse().group("key")

    def parse(self):
        m = re.search("(?P<bucket>[^/]*)/(?P<key>.*)", self.local_file())
        if len(m.groups()) != 2:
            raise WorkflowError(
                "GS remote file {} does not have the form "
                "<bucket>/<key>.".format(self.local_file())
            )
        return m
