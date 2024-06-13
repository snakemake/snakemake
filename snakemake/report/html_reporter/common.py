import os
from pathlib import Path
from snakemake.common import is_local_file
from snakemake.report.common import data_uri_from_file

from snakemake_interface_common.exceptions import WorkflowError


def get_resource_as_string(path_or_uri):
    import requests

    if is_local_file(path_or_uri):
        return open(Path(__file__).parent / "template" / path_or_uri).read()
    else:
        r = requests.get(path_or_uri)
        if r.status_code == requests.codes.ok:
            return r.text
        raise WorkflowError(
            "Failed to download resource needed for report: {}".format(path_or_uri)
        )


def get_result_uri(result, mode_embedded):
    if mode_embedded:
        return data_uri_from_file(result.path)
    else:
        return os.path.join("data/raw", result.id, result.filename)
