__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import email
import os

# module-specific
from snakemake.remote.S3 import RemoteObject, RemoteProvider as S3RemoteProvider


try:
    from google.cloud import storage
except ImportError as e:
    raise WorkflowError("The Python 3 package 'google-cloud-sdk' "
        "needs to be installed to use GS remote() file functionality. %s" % e.msg)


class RemoteProvider(AbstractRemoteProvider):
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

        m = re.search("(?P<bucket>[^/]*)/(?P<key>.*)", self.local_file())
        if len(m.groups()) != 2:
            raise WorkflowError("GS remote file {} does not have the form "
                "<bucket>/<key>.".format(self.local_file()))
        self.key = m.group("key")
        self.bucket_name = m.group("bucket")
        self.bucket = self.client.bucket(self.bucket_name)
        self.blob = self.bucket.blob(self.key)

    # === Implementations of abstract class members ===

    def exists(self):
        return self.blob.exists()

    def mtime(self):
        if self.exists():
            t = self.blob.updated
            # email.utils parsing of timestamp mirrors boto whereas
            # time.strptime() can have TZ issues due to DST
            modified_tuple = email.utils.parsedate_tz(t)
            epoch_time = int(email.utils.mktime_tz(modified_tuple))
            return epoch_time
        else:
            raise WorkflowError("The file does not seem to exist remotely: %s" % self.local_file())

    def size(self):
        if self.exists():
            return self.blob.size // 1024
        else:
            return self._iofile.size_local

    def download(self):
        if self.exists():
            self.blob.download_to_file(self.local_file())
            os.sync()
            return self.local_file()
        return None

    def upload(self):
        self.blob.upload_from_filename(self.local_file())

    @property
    def list(self):
        return [k.name for k in self.bucket.list_blobs()]
