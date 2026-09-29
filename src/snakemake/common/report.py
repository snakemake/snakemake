# compat: python 3.7 (script support)
# Must be kept compatible to Python 3.7 because it is used in Snakemake's
# Python script support.
# Only modules from python standard library or other Snakemake modules that are
# Python 3.7 compatible should be imported here (except for methods
# that are only called by Snakemake itself)!

import hashlib
from pathlib import Path
from typing import Union


def get_report_id(path: Union[str, Path]) -> str:
    h = hashlib.sha256()
    h.update(str(path).encode())

    return h.hexdigest()
