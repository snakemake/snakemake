import os
import sys
import subprocess

sys.path.insert(0, os.path.dirname(__file__))

from common import *


def test_tes():
    subprocess.call(["rm", "-rf", "tests/test_tes/.snakemake"])
    subprocess.call(["rm", "-rf", "tests/test_tes/output.txt"])
    subprocess.call(["rm", "-rf", "tests/test_tes/output.txt.bz2"])
    subprocess.call(["rm", "-rf", "tests/test_tes/test.log"])
    workdir = dpath("test_tes")
    run(
        workdir,
        snakefile="Snakefile",
        tes="http://localhost:8000",
        no_tmpdir=True,
        cleanup=False,
    )
