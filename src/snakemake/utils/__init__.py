# compat: python 3.7 (script support)
# Must be kept compatible to Python 3.7 because it is used in Snakemake's
# Python script support.
# Only modules from python standard library or other Snakemake modules that are
# Python 3.7 compatible should be imported here (except for methods
# that are only called by Snakemake itself)!

# pyrefly: ignore-errors[unused-import]

import sys

# re-export functions from submodules
from snakemake.utils.format import format
from snakemake.utils.argvquote import argvquote, cmd_exe_quote

if sys.version_info >= (3, 11):
    # import utils that are not <3.11 compatible
    from snakemake.utils.paramspace import Paramspace
    from snakemake.utils.validate import validate
    from snakemake.utils.misc import (
        simplify_path,
        linecount,
        listfiles,
        makedirs,
        report,
        R,
        read_job_properties,
        min_version,
        update_config,
        available_cpu_count,
        os_sync,
        find_bash_on_windows,
    )
