import json
import pytest
from pathlib import Path
from .common import run
from snakemake_interface_executor_plugins.settings import DeploymentMethod


# See https://github.com/snakemake/snakemake/pull/3636.
SNAKEFILE = """\
from pathlib import Path

DIR = Path("scripts")

rule all:
    conda: "env.yaml"
    notebook: DIR / "process_data.ipynb"
"""

ENV_YAML = """\
channels:
  - conda-forge
  - bioconda
dependencies:
  - python >=3.5
  - jupyter
  - papermill
"""

NOTEBOOK = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "cells": [],
    "metadata": {
        "kernelspec": {"name": "python3", "display_name": "Python 3"},
        "language_info": {"name": "python"},
    },
}


@pytest.fixture
def testdir(tmpdir):
    p = tmpdir.mkdir("testdir")
    p.mkdir("scripts")
    return p


@pytest.fixture
def testdir_notebook_pathlike(testdir):
    p = testdir.join("scripts", "process_data.ipynb")
    p.write(json.dumps(NOTEBOOK))
    p = testdir.join("env.yaml")
    p.write(ENV_YAML)
    p = testdir.join("Snakefile")
    p.write(SNAKEFILE)
    return testdir


def test_jupyter_notebook_pathlike(testdir_notebook_pathlike):
    run(
        Path(testdir_notebook_pathlike),
        check_results=False,
        deployment_method={DeploymentMethod.CONDA},
    )
