import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from common import run, dpath
from snakemake.common import __version__


def test_issue413():
    run(dpath("test_issue413"), no_tmpdir=True)

    os.unlink(dpath("test_issue413/output.txt"))
    os.unlink(dpath("test_issue413/stats.txt"))
