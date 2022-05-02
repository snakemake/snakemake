import os
import sys
import shlex

sys.path.insert(0, os.path.dirname(__file__))

from common import *

def test_slurm():
    workdir = dpath("test_slurm")
    # first, we need to copy the Snakefile to the "cluster"
    cp_cmd = "vagrant scp Snakefile controller:."
    # this is the command to be tested:
    cmd = "vagrant ssh controller -- -t 'snakemake --use-envmodules --slurm -j 1 -F'"
    # if our workflow was sucessful the 'pi.calc'-file should be present.
    check_cmd = "vagrant scp controller:pi.calc ."

    subprocess.check_call(shlex.split(cp_cmd))
    subprocess.check_call(shlex.split(cmd))
    subprocess.check_call(shlex.split(check_cmd))
