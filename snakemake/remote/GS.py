__author__ = "Johannes Köster"
__copyright__ = "Copyright 2017-2019, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import base64
import os
import re
import struct
import time

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError, CheckSumMismatchException
from snakemake.common import lazy_property
import snakemake.io
from snakemake.utils import os_sync

try:
    import google.cloud
    from google.cloud import storage
    from google.api_core import retry
    from google_crc32c import Checksum
except ImportError as e:
    raise WorkflowError(
        "The Python 3 packages 'google-cloud-sdk' and `google-crc32c` "
        "need to be installed to use GS remote() file functionality. %s" % e.msg
    )


def google_cloud_retry_predicate(ex):
    """Given an exception from Google Cloud, determine if it's one in the
    listing of transient errors (determined by function
    google.api_core.retry.if_transient_error(exception)) or determine if
    triggered by a hash mismatch due to a bad download. This function will
    return a boolean to indicate if retry should be done, and is typically
    used with the google.api_core.retry.Retry as a decorator (predicate).

    Arguments:
      ex (Exception) : the exception passed from the decorated function
    Returns: boolean to indicate doing retry (True) or not (False)
    """
    # Most likely case is Google API transient error
    if retry.if_transient_error(ex):
        return True
    # Could also be checksum mismatch of download
    if isinstance(ex, CheckSumMismatchException):
        return True
    return False


@retry.Retry(predicate=google_cloud_retry_predicate)
def download_blob(blob, filename):
    """A helper function to download a storage Blob to a blob_file (the filename)
    and validate it using the Crc32cCalculator.

    Arguments:
      blob (storage.Blob) : the Google storage blob object
      blob_file (str)     : the file path to download to
    Returns: boolean to indicate doing retry (True) or not (False)
    """

    # create parent directories if necessary
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    # ideally we could calculate hash while streaming to file with provided function
    # https://github.com/googleapis/python-storage/issues/29
    with open(filename, "wb") as blob_file:
        parser = Crc32cCalculator(blob_file)
        blob.download_to_file(parser)
    os.sync()

    # **Important** hash can be incorrect or missing if not refreshed
    blob.reload()

    # Compute local hash and verify correct
    if parser.hexdigest() != blob.crc32c:
        os.remove(filename)
        raise CheckSumMismatchException("The checksum of %s does not match." % filename)
    return filename


class Crc32cCalculator:
    """The Google Python client doesn't provide a way to stream a file being
    written, so we can wrap the file object in an additional class to
    do custom handling. This is so we don't need to download the file
    and then stream read it again to calculate the hash.
    """

    def __init__(self, fileobj):
        self._fileobj = fileobj
        self.checksum = Checksum()

    def write(self, chunk):
        self._fileobj.write(chunk)
        self._update(chunk)

    def _update(self, chunk):
        """Given a chunk from the read in file, update the hexdigest"""
        self.checksum.update(chunk)

    def hexdigest(self):
        """Return the hexdigest of the hasher.
        The Base64 encoded CRC32c is in big-endian byte order.
        See https://cloud.google.com/storage/docs/hashes-etags
        """
        return base64.b64encode(self.checksum.digest()).decode("utf-8")


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

    async def inventory(self, cache: snakemake.io.IOCache):
        """Using client.list_blobs(), we want to iterate over the objects in
        the "folder" of a bucket and store information about the IOFiles in the
        provided cache (snakemake.io.IOCache) indexed by bucket/blob name.
        This will be called by the first mention of a remote object, and
        iterate over the entire bucket once (and then not need to again).
        This includes:
         - cache.exist_remote
         - cache_mtime
         - cache.size
        """
        if cache.remaining_wait_time <= 0:
            # No more time to create inventory.
            return

        start_time = time.time()
        subfolder = os.path.dirname(self.blob.name)
        for blob in self.client.list_blobs(self.bucket_name, prefix=subfolder):
            # By way of being listed, it exists. mtime is a datetime object
            name = "{}/{}".format(blob.bucket.name, blob.name)
            cache.exists_remote[name] = True
            cache.mtime[name] = blob.updated.timestamp()
            cache.size[name] = blob.size

        cache.remaining_wait_time -= time.time() - start_time

        # Mark bucket and prefix as having an inventory, such that this method is
        # only called once for the subfolder in the bucket.
        cache.exists_remote.has_inventory.add("%s/%s" % (self.bucket_name, subfolder))

    # === Implementations of abstract class members ===

    def get_inventory_parent(self):
        return self.bucket_name

    @retry.Retry(predicate=google_cloud_retry_predicate)
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

    @retry.Retry(predicate=google_cloud_retry_predicate, deadline=600)
    def download(self):
        """Download with maximum retry duration of 600 seconds (10 minutes)"""
        if not self.exists():
            return None

        # Create just a directory, or a file itself
        if snakemake.io.is_flagged(self.local_file(), "directory"):
            return self._download_directory()
        return download_blob(self.blob, self.local_file())

    @retry.Retry(predicate=google_cloud_retry_predicate)
    def _download_directory(self):
        """A 'private' function to handle download of a storage folder, which
        includes the content found inside.
        """
        # Create the directory locally
        os.makedirs(self.local_file(), exist_ok=True)

        for blob in self.client.list_blobs(self.bucket_name, prefix=self.key):
            local_name = "{}/{}".format(blob.bucket.name, blob.name)

            # Don't try to create "directory blob"
            if os.path.exists(local_name) and os.path.isdir(local_name):
                continue

            download_blob(blob, local_name)

        # Return the root directory
        return self.local_file()

    @retry.Retry(predicate=google_cloud_retry_predicate)
    def upload(self):
        try:
            if not self.bucket.exists():
                self.bucket.create()
                self.update_blob()

            # Distinguish between single file, and folder
            f = self.local_file()
            if os.path.isdir(f):

                # Ensure the "directory" exists
                self.blob.upload_from_string(
                    "", content_type="application/x-www-form-urlencoded;charset=UTF-8"
                )
                for root, _, files in os.walk(f):
                    for filename in files:
                        filename = os.path.join(root, filename)
                        bucket_path = filename.lstrip(self.bucket.name).lstrip("/")
                        blob = self.bucket.blob(bucket_path)
                        blob.upload_from_filename(filename)
            else:
                self.blob.upload_from_filename(f)
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

    @retry.Retry(predicate=google_cloud_retry_predicate)
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

    @property
    def key(self):
        key = self.parse().group("key")
        f = self.local_file()
        if snakemake.io.is_flagged(f, "directory"):
            key = key if f.endswith("/") else key + "/"
        return key

    def parse(self):
        m = re.search("(?P<bucket>[^/]*)/(?P<key>.*)", self.local_file())
        if len(m.groups()) != 2:
            raise WorkflowError(
                "GS remote file {} does not have the form "
                "<bucket>/<key>.".format(self.local_file())
            )
        return m
