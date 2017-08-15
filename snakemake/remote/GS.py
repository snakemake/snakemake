__author__ = "Johannes Köster"
__copyright__ = "Copyright 2017, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import os
import re

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError
from snakemake.common import lazy_property

try:
    import google.cloud
    from google.cloud import storage
except ImportError as e:
    raise WorkflowError("The Python 3 package 'google-cloud-sdk' "
        "needs to be installed to use GS remote() file functionality. %s" % e.msg)


class RemoteProvider(AbstractRemoteProvider):

    supports_default = True

    def __init__(self, *args, stay_on_remote=False, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)

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
    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)

        if provider:
            self.client = provider.remote_interface()
        else:
            self.client = storage.Client(*args, **kwargs)

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
            raise WorkflowError("The file does not seem to exist remotely: %s" % self.local_file())

    def size(self):
        if self.exists():
            self.update_blob()
            return self.blob.size // 1024
        else:
            return self._iofile.size_local

    def download(self):
        if self.exists():
            os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)
            self.blob.download_to_filename(self.local_file())
            os.sync()
            return self.local_file()
        return None

    def upload(self):
        try:
            if not self.bucket.exists():
                self.bucket.create()
                self.update_blob()
            self.blob.upload_from_filename(self.local_file())
        except google.cloud.exceptions.Forbidden as e:
            raise WorkflowError(e,
                "When running locally, make sure that you are authenticated "
                "via gcloud (see Snakemake documentation). When running in a "
                "kubernetes cluster, make sure that storage-rw is added to "
                "--scopes (see Snakemake documentation).")

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
        return self.client.bucket(self.bucket_name)

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
            raise WorkflowError("GS remote file {} does not have the form "
                "<bucket>/<key>.".format(self.local_file()))
        return m
