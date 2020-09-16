import os
import sys
import tempfile

from google.cloud import storage

sys.path.insert(0, os.path.dirname(__file__))

from common import *


def has_google_credentials():
    return "GOOGLE_APPLICATION_CREDENTIALS" in os.environ


google_credentials = pytest.mark.skipif(
    not has_google_credentials(),
    reason="Skipping google lifesciences tests because  "
    "GOOGLE_APPLICATION_CREDENTIALS not found in the environment.",
)


def cleanup_google_storage(prefix, bucket_name="snakemake-testing"):
    """Given a storage prefix and a bucket, recursively delete files there
    This is intended to run at the completion of testing to ensure that
    the bucket is cleaned up.

    Arguments:
      prefix (str) : the "subfolder" or prefix for some files in the buckets
      bucket_name (str) : the name of the bucket, default snakemake-testing
    """
    client = storage.Client()
    bucket = client.get_bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=prefix)
    for blob in blobs:
        blob.delete()


@google_credentials
def test_google_lifesciences():
    storage_prefix = "snakemake-testing-%s" % next(tempfile._get_candidate_names())
    workdir = dpath("test_google_lifesciences")
    try:
        run(
            workdir,
            use_conda=True,
            default_remote_prefix="snakemake-testing/%s" % storage_prefix,
            google_lifesciences=True,
            google_lifesciences_cache=True,
            preemption_default=None,
            preemptible_rules=["pack=1"],
        )
    finally:
        cleanup_google_storage(storage_prefix)
