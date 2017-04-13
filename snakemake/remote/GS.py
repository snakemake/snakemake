__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# module-specific
from snakemake.remote.S3 import RemoteObject, RemoteProvider as S3RemoteProvider


class RemoteProvider(S3RemoteProvider):
    def __init__(self, *args, stay_on_remote=False, **kwargs):
        kwargs["host"] = "storage.googleapis.com"
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return 'gs://'

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ['s3://', 'gs://']
