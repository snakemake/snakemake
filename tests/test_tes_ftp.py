import os
import sys
import subprocess

sys.path.insert(0, os.path.dirname(__file__))

from common import *

def test_tes_ftp():
    subprocess.call(["rm", "-rf", "tests/test_tes_ftp/.snakemake"])
    subprocess.call(["rm", "-rf", "tests/test_tes_ftp/output.txt"])
    subprocess.call(["rm", "-rf", "tests/test_tes_ftp/output.txt.bz2"])
    subprocess.call(["rm", "-rf", "tests/test_tes_ftp/stats.txt"])
    subprocess.call(["rm", "-rf", "tests/test_tes_ftp/test.log"])
    workdir = dpath("test_tes_ftp")
    run(
        workdir,
        snakefile="Snakefile",
        tes="http://localhost:8000",
        use_conda=True,
        conda_prefix="/tmp/conda",
        conda_frontend="conda",
        envvars=["HTTP_PROXY", "HTTPS_PROXY", "CONDA_PKGS_DIRS", "CONDA_ENVS_PATH"],
        no_tmpdir=True,
        cleanup=False
    )
