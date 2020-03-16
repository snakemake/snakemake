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

    envfile = lint.joinpath(".env.py")
    print(envfile)
    if envfile.exists():
        exec(open(envfile).read())
        print(envfile)
    out = (
        sp.check_output(
            [
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
        assert out == ""
    else:
        print(out)
        assert out
