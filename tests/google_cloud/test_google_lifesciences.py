from . import *


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
        )
    finally:
        cleanup_google_storage(storage_prefix, bucket_name)
