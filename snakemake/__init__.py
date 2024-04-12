__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

import sys
from snakemake.common import __version__, PIP_DEPLOYMENTS_PATH

sys.path.append(PIP_DEPLOYMENTS_PATH)

# Reexports that are part of the public API:
from snakemake.shell import shell


if __name__ == "__main__":
    from snakemake.cli import main
    import sys

    main(sys.argv)
