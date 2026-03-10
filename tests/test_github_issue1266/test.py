#!/usr/bin/env python3

import sys

# In principle,
# the fact that the script is able to run at all
# means that issue #1266 is fixed,
# since in Snakemake 6.10 the crash occurs
# while running the setup code that Snakemake injects.
# However,
# for good measure,
# we also ensure that we have access to
# the packages associated with the Snakemake-managed conda environment,
# so we are confident that the environment is available.
assert any(
    ".snakemake/conda" in path.replace("\\", "/")
    for path in sys.path
    if path.endswith("site-packages")
)
