import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path

def test_b():
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        shutil.copytree(".tests/unit/b", workdir)
        sp.check_output([
            "snakemake", 
            "test/0.tsv", 
            "-F", 
            "-j1",
            "--directory",
            workdir,
        ])
