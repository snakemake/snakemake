import sys
import os
import re

sys.path.insert(0, os.path.dirname(__file__))

from common import *


def test_az_batch_executor():
    # AZ_BATCH_ACCOUNT_URL=https://${batch_account_name}.${region}.batch.azure.com
    bau = os.getenv("AZ_BATCH_ACCOUNT_URL")
    prefix = os.getenv("AZ_BLOB_PREFIX")
    wdir = dpath("test_azure_batch")
    blob_account_url = os.getenv("AZ_BLOB_ACCOUNT_URL")
    assert blob_account_url is not None and blob_account_url.strip() != ""

    sas_token = os.getenv("AZ_BLOB_CREDENTIAL") 
    assert sas_token is not None, "SAS_TOKEN environment variable is not set."

    # pattern for Storage Account SAS Token
    pattern_account = re.compile(r"\?sv=.*&ss=.*&srt=.*&sp=.*&se=.*&st=.*&spr=.*&sig=.*")

    assert pattern_account.match(sas_token) is not None, "AZ_BLOB_CREDENTIAL does not match the SAS token pattern."

    run(
        path=wdir,
        default_remote_prefix=prefix,
        container_image="jakevc/snakemake",
        envvars=["AZ_BLOB_ACCOUNT_URL", "AZ_BLOB_CREDENTIAL"],
        az_batch=True,
        az_batch_account_url=bau,
    )
