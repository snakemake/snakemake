# compat: python 3.7 (script support)
# Must be kept compatible to Python 3.7 because it is used in Snakemake's
# Python script support.
# Only modules from python standard library or other Snakemake modules that are
# Python 3.7 compatible should be imported here (except for methods
# that are only called by Snakemake itself)!

from pathlib import Path
import platform
import uuid
from typing import Tuple

RULEFUNC_CONTEXT_MARKER = "__is_snakemake_rule_func"
MIN_PY_VERSION: Tuple[int, int] = (3, 7)
UUID_NAMESPACE = uuid.uuid5(uuid.NAMESPACE_URL, "https://snakemake.readthedocs.io")
NOTHING_TO_BE_DONE_MSG = (
    "Nothing to be done (all requested files are present and up to date)."
)
ON_WINDOWS = platform.system() == "Windows"
# limit the number of input/output files list in job properties
# see https://github.com/snakemake/snakemake/issues/2097
IO_PROP_LIMIT = 100
SNAKEFILE_CHOICES = list(
    map(
        Path,
        (
            "Snakefile",
            "snakefile",
            "workflow/Snakefile",
            "workflow/snakefile",
        ),
    )
)
