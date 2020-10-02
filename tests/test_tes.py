import os
import sys
import subprocess

sys.path.insert(0, os.path.dirname(__file__))

from common import *


def test_tes():
    workdir = dpath("test_tes")
    run(
        workdir,
        snakefile="Snakefile",
        tes="http://localhost:8000",
        no_tmpdir=True,
        cleanup=False,
    )
