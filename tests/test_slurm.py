import os
import sys
import shlex

from common import *

def test_slurm():
    run(dpath("test_slurm"), slurm=True)
