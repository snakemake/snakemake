import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from common import *


def test_az_batch_executor():
    bau = os.getenv("AZ_BATCH_ACCOUNT_URL")
    wdir = dpath("test_azure_batch")
    run(
        path=wdir,
        default_remote_prefix="snaketest",
        container_image="jakevc/snakemake",
        envvars=["AZ_BLOB_ACCOUNT_URL", "AZ_BLOB_CREDENTIAL"],
        az_batch=True,
        az_batch_account_url=bau,
    )
