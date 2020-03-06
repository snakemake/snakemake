import os
import sys
import subprocess

sys.path.insert(0, os.path.dirname(__file__))

from common import *

# run tibanna test before any moto-related tests because they apparently render AWS environment variables invalid or uneffective.
def test_tes():
    subprocess.call(["rm", "-rf", "tests/test_tes/.snakemake"])
    subprocess.call(["rm", "-rf", "tests/test_tes/output.txt"])
    workdir = dpath("test_tes")
    run(
        workdir,
        tes=True,
        no_tmpdir=True,
        cleanup=False
    )
