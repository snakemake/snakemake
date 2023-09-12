import os
import requests
import sys
import tempfile
import google.auth

from google.cloud import storage

sys.path.insert(0, os.path.dirname(__file__))

from common import *


def has_google_credentials():
    credentials, _ = google.auth.default()
    return credentials


google_credentials = pytest.mark.skipif(
    not has_google_credentials(),
    reason="Skipping google lifesciences tests because  "
    "Google credentials were not found in the environment.",
)


def get_default_service_account_email():
    """Returns the default service account if running on a GCE VM, otherwise None."""
    response = requests.get(
        "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/email",
        headers={"Metadata-Flavor": "Google"},
    )
    if response.status_code == requests.codes.ok:
        return response.text
    else:
        return None


def cleanup_google_storage(prefix, bucket_name="snakemake-testing", restrict_to=None):
    """Given a storage prefix and a bucket, recursively delete files there
    This is intended to run after testing to ensure that
    the bucket is cleaned up.

    Arguments:
      prefix (str) : the "subfolder" or prefix for some files in the buckets
      bucket_name (str) : the name of the bucket, default snakemake-testing
      restrict_to (list) : only delete files in these paths (None deletes all)
    """
    client = storage.Client()
    bucket = client.get_bucket(bucket_name)
    blobs = bucket.list_blobs(prefix="source")
    for blob in blobs:
        blob.delete()
    blobs = bucket.list_blobs(prefix=prefix)
    for blob in blobs:
        if restrict_to is None or f"{bucket_name}/{blob.name}" in restrict_to:
            blob.delete()
    if restrict_to is None:
        # Using API we get an exception about bucket deletion
        shell("gsutil -m rm -r gs://{bucket.name}/* || true")
        bucket.delete()


def create_google_storage(bucket_name="snakemake-testing"):
    """Given a bucket name, create the Google storage bucket,
    intending to be used for testing and then cleaned up by
    cleanup_google_storage

    Arguments:
      bucket_name (str) : the name of the bucket, default snakemake-testing
    """
    client = storage.Client()
    return client.create_bucket(bucket_name)


def get_temp_bucket_name():
    return "snakemake-testing-%s-bucket" % next(tempfile._get_candidate_names())


@google_credentials
def test_google_lifesciences():
    bucket_name = get_temp_bucket_name()
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
            google_lifesciences_service_account_email=get_default_service_account_email(),
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
    bucket_name = get_temp_bucket_name()
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
            google_lifesciences_service_account_email=get_default_service_account_email(),
        )
    finally:
        cleanup_google_storage(storage_prefix, bucket_name)


@google_credentials
def test_cloud_checkpoints_issue574():
    """see Github issue #574"""
    bucket_name = get_temp_bucket_name()
    create_google_storage(bucket_name)
    storage_prefix = "test_cloud_checkpoints_issue574"
    workdir = dpath("test_cloud_checkpoints_issue574")
    try:
        run(
            workdir,
            use_conda=True,
            default_remote_prefix="%s/%s" % (bucket_name, storage_prefix),
            google_lifesciences=True,
            google_lifesciences_cache=False,
            google_lifesciences_service_account_email=get_default_service_account_email(),
        )
    finally:
        cleanup_google_storage(storage_prefix, bucket_name)


def test_github_issue1396():
    bucket_name = get_temp_bucket_name()
    create_google_storage(bucket_name)
    storage_prefix = "test_github_issue1396"
    workdir = dpath("test_github_issue1396")
    try:
        run(
            workdir,
            default_remote_prefix="%s/%s" % (bucket_name, storage_prefix),
            google_lifesciences=True,
            google_lifesciences_cache=False,
            dryrun=True,
        )
    finally:
        cleanup_google_storage(storage_prefix, bucket_name)


def test_github_issue1460():
    service_account_email = get_default_service_account_email()
    bucket_name = get_temp_bucket_name()
    create_google_storage(bucket_name)
    storage_prefix = "test_github_issue1460"
    prefix = "%s/%s" % (bucket_name, storage_prefix)
    workdir = dpath("test_github_issue1460")
    try:
        run(
            workdir,
            default_remote_prefix=prefix,
            google_lifesciences=True,
            google_lifesciences_cache=False,
            google_lifesciences_service_account_email=service_account_email,
        )
        cleanup_google_storage(
            storage_prefix,
            bucket_name,
            restrict_to=[
                f"{prefix}/test.txt",
                f"{prefix}/blob.txt",
                f"{prefix}/pretest.txt",
            ],
        )
        run(
            workdir,
            default_remote_prefix=prefix,
            google_lifesciences=True,
            google_lifesciences_cache=False,
            google_lifesciences_service_account_email=service_account_email,
        )
    finally:
        cleanup_google_storage(storage_prefix, bucket_name)
