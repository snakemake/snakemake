import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path

def test_a():
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        shutil.copytree(".tests/unit/a", workdir)
        sp.check_output([
            "snakemake", 
            "test/0.txt", 
            "-F", 
            "-j1",
            "--directory",
            workdir,
        ])
