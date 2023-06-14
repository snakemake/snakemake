import sys
import os
import re

sys.path.insert(0, os.path.dirname(__file__))

from common import *


@azbatch
def test_az_batch_executor():
    # AZ_BATCH_ACCOUNT_URL=https://${batch_account_name}.${region}.batch.azure.com
    bau = os.getenv("AZ_BATCH_ACCOUNT_URL")
    prefix = os.getenv("AZ_BLOB_PREFIX")
    wdir = dpath("test_azure_batch")
    blob_account_url = os.getenv("AZ_BLOB_ACCOUNT_URL")
    assert blob_account_url is not None and blob_account_url.strip() != ""

    run(
        path=wdir,
        default_remote_prefix=prefix,
        container_image="snakemake/snakemake",
        envvars=["AZ_BLOB_ACCOUNT_URL", "AZ_BLOB_CREDENTIAL"],
        az_batch=True,
        az_batch_account_url=bau,
    )
