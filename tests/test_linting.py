import os
from pathlib import Path
import subprocess as sp
from itertools import product

import pytest

LINT_DIR = Path(__file__).parent.joinpath("linting")


@pytest.mark.parametrize(
    "lint, case", product(os.listdir(LINT_DIR), ["positive", "negative"])
)
def test_lint(lint, case):
    lint = LINT_DIR.joinpath(lint)

    try:
        out = (
            sp.check_output(
                [
                    "python",
                    "-m",
                    "snakemake",
                    "--lint",
                    "--directory",
                    lint,
                    "--snakefile",
                    lint.joinpath(case).with_suffix(".smk"),
                ],
                stderr=sp.STDOUT,
            )
            .decode()
            .strip()
        )
        if case == "positive":
            assert out == "Congratulations, your workflow is in a good condition!"
        else:
            print(out)
            assert (
                False
            ), "Negative lint example but linting command exited with status 0."

    except sp.CalledProcessError as e:
        if case == "negative":
            assert e.output.decode().strip()
        else:
            raise e

    else:
        print(out)
        assert out
