import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path

def test_all():
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        shutil.copytree(".tests/unit/all", workdir)
        sp.check_output([
            "snakemake", 
            "all", 
            "-F", 
            "-j1",
            "--directory",
            workdir,
        ])
