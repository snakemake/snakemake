__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# module-specific
from snakemake.remote.S3 import RemoteObject, RemoteProvider as S3RemoteProvider

try:
    # third-party modules
    import boto
    from boto.s3.key import Key
    from filechunkio import FileChunkIO
except ImportError as e:
    raise WorkflowError("The Python 3 packages 'boto' and 'filechunkio' " + 
        "need to be installed to use S3 remote() file functionality. %s" % e.msg)

class RemoteProvider(S3RemoteProvider):
    def __init__(self, *args, **kwargs):
        if "gs_access_key_id" in kwargs:
            kwargs["aws_access_key_id"] = kwargs.pop("gs_access_key_id")
        if "gs_secret_access_key" in kwargs:
            kwargs["aws_secret_access_key"] = kwargs.pop("gs_secret_access_key")

        kwargs["host"] = "storage.googleapis.com"
        super(RemoteProvider, self).__init__(*args, **kwargs)
