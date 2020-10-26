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
    blobs = bucket.list_blobs(prefix="source")
    for blob in blobs:
        blob.delete()
    # Using API we get an exception about bucket deletion
    shell("gsutil -m rm -r gs://{bucket.name}/* || true")
    bucket.delete()


def create_google_storage(bucket_name="snakemake-testing"):
    """Given a bucket name, create the Google storage bucket,
    with intention to be used for testing and then cleaned up by
    cleanup_google_storage

    Arguments:
      bucket_name (str) : the name of the bucket, default snakemake-testing
    """
    client = storage.Client()
    return client.create_bucket(bucket_name)


@google_credentials
def test_google_lifesciences():
    bucket_name = "snakemake-testing-%s" % next(tempfile._get_candidate_names())
    create_google_storage(bucket_name)
    storage_prefix = "test_google_lifesciences"
    workdir = dpath("test_google_lifesciences")
    try:
        run(
            workdir,
            use_conda=True,
            default_remote_prefix="%s/%s" % (bucket_name, storage_prefix),
            google_lifesciences=True,
            google_lifesciences_cache=False,
            preemption_default=None,
            preemptible_rules=["pack=1"],
        )
    finally:
        cleanup_google_storage(storage_prefix, bucket_name)


@pytest.mark.skip(
    reason="Cannot test using touch with a remote prefix until the container image is deployed."
)
@google_credentials
def test_touch_remote_prefix():
    bucket_name = "snakemake-testing-%s" % next(tempfile._get_candidate_names())
    create_google_storage(bucket_name)
    storage_prefix = "test_touch_remote_prefix"
    workdir = dpath("test_touch_remote_prefix")
    try:
        run(
            workdir,
            use_conda=True,
            default_remote_prefix="%s/%s" % (bucket_name, storage_prefix),
            google_lifesciences=True,
            google_lifesciences_cache=False,
        )
    finally:
        cleanup_google_storage(storage_prefix, bucket_name)
