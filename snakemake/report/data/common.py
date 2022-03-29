from pathlib import Path
from snakemake.common import is_local_file
from snakemake.exceptions import WorkflowError


def get_resource_as_string(path_or_uri):
    import requests

    if is_local_file(path_or_uri):
        return open(Path(__file__).parent.parent / "template" / path_or_uri).read()
    else:
        r = requests.get(path_or_uri)
        if r.status_code == requests.codes.ok:
            return r.text
        raise WorkflowError(
            "Failed to download resource needed for " "report: {}".format(path_or_uri)
        )
