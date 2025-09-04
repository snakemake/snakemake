"""
Common code for unit testing of rules generated with Snakemake 9.8.2.dev50.
"""

import os
from pathlib import Path
from subprocess import check_output


class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir

    def check(self):
        # Input files
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        print(f"input: {input_files}")  # DEBUG
        # Workdir files
        workdir_files = set(
            (Path(path) / f).relative_to(self.workdir)
            for path, subdirs, files in os.walk(self.workdir)
            for f in files
        )
        print(f"workdir: {workdir_files}")  # DEBUG
        # Expected files
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        print(f"expected: {expected_files}")  # DEBUG

        assert expected_files.issubset(
            workdir_files
        ), f"Output files missing: {expected_files - workdir_files}"

        # Compare output and expected files
        for f in expected_files:
            self.compare_files(self.expected_path / f, self.workdir / f)

    def compare_files(self, expected_file, generated_file):
        if expected_file.suffix == ".gz":
            cmp = "zcmp"
        elif expected_file.suffix == ".bz2":
            cmp = "bzcmp"
        elif expected_file.suffix == ".xz":
            cmp = "xzcmp"
        else:
            cmp = "cmp"

        check_output([cmp, expected_file, generated_file])
