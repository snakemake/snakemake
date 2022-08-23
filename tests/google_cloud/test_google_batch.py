from . import *


@google_credentials
def test_google_batch():
    bucket_name = get_temp_bucket_name()
    create_google_storage(bucket_name)
    storage_prefix = "test_google_batch"
    workdir = dpath("test_google_lifesciences")
    try:
        run(
            workdir,
            use_conda=True,
            default_remote_prefix=f"{bucket_name}/{storage_prefix}",
            google_batch=True,
            preemption_default=None,
            preemptible_rules=["pack=1"],
        )
    finally:
        cleanup_google_storage(storage_prefix, bucket_name)
