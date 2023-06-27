import os
import boto3
import base64
import uuid
import requests


def detect_aws_region():
    # check environment variables
    for ev in ("AWS_REGION", "AWS_REGION"):
        if os.environ.get(ev):
            return os.environ[ev]

    # check boto3, which will load ~/.aws
    if boto3.DEFAULT_SESSION and boto3.DEFAULT_SESSION.region_name:
        return boto3.DEFAULT_SESSION.region_name
    session = boto3.Session()
    if session.region_name:
        return session.region_name

    # query EC2 metadata
    try:
        return requests.get(
            "http://169.254.169.254/latest/meta-data/placement/region", timeout=2.0
        ).text
    except:
        pass

    return None


def efs_id_from_access_point(region_name, fsap_id):
    # Resolve the EFS access point id (fsap-xxxx) to the associated file system id (fs-xxxx). Saves
    # user from having to specify both.
    aws_efs = boto3.Session().client("efs", region_name=region_name)
    desc = aws_efs.describe_access_points(AccessPointId=fsap_id)
    assert len(desc.get("AccessPoints", [])) == 1
    desc = desc["AccessPoints"][0]
    fs_id = desc["FileSystemId"]
    assert isinstance(fs_id, str) and fs_id.startswith("fs-")
    return fs_id
