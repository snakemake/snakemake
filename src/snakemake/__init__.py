__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

import sys

# Reexports that are part of the public API:
from snakemake._version import version as __version__
from snakemake.shell import shell

PIP_DEPLOYMENTS_PATH = ".snakemake/pip-deployments"

sys.path.append(PIP_DEPLOYMENTS_PATH)

if __name__ == "__main__":
    import sys

    from snakemake.cli import main

    main(sys.argv)
