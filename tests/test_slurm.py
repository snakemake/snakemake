__authors__ = ["Christian Meesters", "Johannes Köster"]
__copyright__ = "Copyright 2022, Christian Meesters, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from .common import *
from .conftest import skip_on_windows


@skip_on_windows
def test_slurm_mpi():
    run(
        dpath("test_slurm_mpi"),
        slurm=True,
        show_failed_logs=True,
        use_conda=True,
        default_resources=DefaultResources(["account=runner", "partition=debug"]),
    )


@skip_on_windows
def test_slurm_group_job():
    """
    same test as test_group_job(),
    but for SLURM - checks whether
    the group-property is correctly
    propagated.
    """
    run(
        dpath("test_group_job"),
        slurm=True,
        default_resources=DefaultResources(["account=runner", "partition=debug"]),
    )


@skip_on_windows
def test_slurm_complex():
    os.environ["TESTVAR"] = "test"
    os.environ["TESTVAR2"] = "test"
    run(
        dpath("test14"),
        snakefile="Snakefile.nonstandard",
        slurm=True,
        default_resources=DefaultResources(["account=runner", "partition=debug"]),
    )