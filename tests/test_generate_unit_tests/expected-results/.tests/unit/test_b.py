"""
Rule test code for unit testing of rules generated with Snakemake 9.8.2.dev2.
"""


import os
import sys
import shutil
import tempfile
import subprocess as sp
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))


def test_b(conda_prefix):

    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = Path(".tests/unit/b/data")
        expected_path = Path(".tests/unit/b/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Copy config/ (if exists) to the temporary workdir.
        config = Path("config")
        if config.exists() and config.is_dir():
            shutil.copytree("config", workdir / config, dirs_exist_ok=True)

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "test/0.tsv",
                "--snakefile",
                "Snakefile",
                "-f",
                "--notemp",
                "-j1",
                "--target-files-omit-workdir-adjustment",
                "--configfile",
                "config.json",
                "--directory",
                workdir,
            ]
            + conda_prefix
        )

        # Check the output byte by byte using cmp/zmp/bzcmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        import common
        common.OutputChecker(data_path, expected_path, workdir).check()
