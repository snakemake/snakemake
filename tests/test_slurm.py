import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from common import *

def test_slurm():
    workdir = dpath("test_slurm")

    run(
        workdir,
        use_conda=False,
        # not using env_modules in this context, either
        # MPI is supposted to be installed, nothing
        # else is needed in the vagrant setup.
        slurm=True,
    )
