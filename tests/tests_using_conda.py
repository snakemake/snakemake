__authors__ = ["Tobias Marschall", "Marcel Martin", "Johannes Köster"]
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import shutil
import sys
import subprocess as sp
from pathlib import Path
import tempfile

import pytest
from snakemake.deployment.conda import get_env_setup_done_flag_file

sys.path.insert(0, os.path.dirname(__file__))

from .common import run, dpath, apptainer, conda, connected
from .conftest import skip_on_windows, only_on_windows, ON_WINDOWS

from snakemake_interface_executor_plugins.settings import (
    DeploymentMethod,
)


xfail_permissionerror_on_win = (
    pytest.mark.xfail(raises=PermissionError) if ON_WINDOWS else lambda x: x
)


@skip_on_windows
@conda
def test_script():
    run(
        dpath("test_script"),
        deployment_method={DeploymentMethod.CONDA},
        check_md5=False,
    )


@skip_on_windows
@conda
def test_script_rs():
    run(
        dpath("test_script_rs"),
        deployment_method={DeploymentMethod.CONDA},
        check_md5=False,
    )


@conda
def test_conda():
    run(dpath("test_conda"), deployment_method={DeploymentMethod.CONDA})


@conda
def test_conda_list_envs():
    run(dpath("test_conda"), conda_list_envs=True, check_results=False)


# TODO failing with FAILED tests/tests.py::test_conda_create_envs_only -
# PermissionError: [WinError 32] The process cannot access the file because
# it is being used by another process:
# 'C:\\Users\\RUNNER~1\\AppData\\Local\\Temp\\snakemake-2q4osog0\\test-env.yaml'
@skip_on_windows
@conda
def test_conda_create_envs_only():
    tmpdir = run(
        dpath("test_conda"),
        conda_create_envs=True,
        check_results=False,
        cleanup=False,
        cleanup_scripts=False,
    )
    env_dir = next(
        (p for p in Path(tmpdir, ".snakemake", "conda").iterdir() if p.is_dir()), None
    )
    assert env_dir is not None
    assert get_env_setup_done_flag_file(Path(env_dir)).exists()
    shutil.rmtree(tmpdir)


@conda
def test_upstream_conda():
    run(
        dpath("test_conda"),
        deployment_method={DeploymentMethod.CONDA},
        conda_frontend="conda",
    )


@skip_on_windows
@conda
def test_deploy_script():
    run(dpath("test_deploy_script"), deployment_method={DeploymentMethod.CONDA})


@skip_on_windows
@conda
def test_deploy_hashing():
    tmpdir = run(
        dpath("test_deploy_hashing"),
        deployment_method={DeploymentMethod.CONDA},
        cleanup=False,
    )
    assert len(next(os.walk(os.path.join(tmpdir, ".snakemake/conda")))[1]) == 2


@conda
def test_conda_custom_prefix():
    run(
        dpath("test_conda_custom_prefix"),
        deployment_method={DeploymentMethod.CONDA},
        conda_prefix="custom",
        set_pythonpath=False,
    )


@only_on_windows
@conda
def test_conda_cmd_exe():
    # Tests the conda environment activation when cmd.exe
    # is used as the shell
    run(dpath("test_conda_cmd_exe"), deployment_method={DeploymentMethod.CONDA})


@skip_on_windows  # wrappers are for linux and macos only
@conda
def test_wrapper():
    run(
        dpath("test_wrapper"),
        deployment_method={DeploymentMethod.CONDA},
        check_md5=False,
    )


@skip_on_windows  # wrappers are for linux and macos only
@conda
def test_wrapper_local_git_prefix():
    import git

    with tempfile.TemporaryDirectory() as tmpdir:
        print("Cloning wrapper repo...")
        repo = git.Repo.clone_from(
            "https://github.com/snakemake/snakemake-wrappers", tmpdir
        )
        print(f"Cloning of {repo} complete.")

        run(
            dpath("test_wrapper"),
            deployment_method={DeploymentMethod.CONDA},
            wrapper_prefix=f"git+file://{tmpdir}",
            check_md5=False,
        )


@skip_on_windows
@apptainer
@connected
@conda
def test_singularity_conda():
    run(
        dpath("test_singularity_conda"),
        deployment_method={DeploymentMethod.CONDA, DeploymentMethod.APPTAINER},
        conda_frontend="conda",
    )


@conda
@pytest.mark.needs_envmodules
def test_archive():
    run(dpath("test_archive"), archive="workflow-archive.tar.gz")


@skip_on_windows
@conda
def test_issue635():
    run(
        dpath("test_issue635"),
        deployment_method={DeploymentMethod.CONDA},
        check_md5=False,
    )


@conda
def test_issue1046():
    run(dpath("test_issue1046"))


@skip_on_windows
@conda
def test_issue1093():
    run(dpath("test_issue1093"), deployment_method={DeploymentMethod.CONDA})


@conda
def test_jupyter_notebook():
    run(dpath("test_jupyter_notebook"), deployment_method={DeploymentMethod.CONDA})


@conda
def test_jupyter_notebook_draft():
    from snakemake.settings.types import NotebookEditMode

    run(
        dpath("test_jupyter_notebook_draft"),
        deployment_method={DeploymentMethod.CONDA},
        edit_notebook=NotebookEditMode(draft_only=True),
        targets=["results/result_intermediate.txt"],
        check_md5=False,
    )


@skip_on_windows
@apptainer
@conda
def test_containerized():
    run(
        dpath("test_containerized"),
        deployment_method={DeploymentMethod.CONDA, DeploymentMethod.APPTAINER},
    )


@skip_on_windows
@conda
def test_containerize():
    run(dpath("test_conda"), containerize=True, check_results=False)


@conda
def test_converting_path_for_r_script():
    run(
        dpath("test_converting_path_for_r_script"),
        cores=1,
        deployment_method={DeploymentMethod.CONDA},
    )


@skip_on_windows
@conda
def test_conda_named():
    run(dpath("test_conda_named"), deployment_method={DeploymentMethod.CONDA})


@skip_on_windows
@conda
def test_conda_function():
    run(
        dpath("test_conda_function"),
        deployment_method={DeploymentMethod.CONDA},
        cores=1,
    )


@skip_on_windows  # the testcase only has a linux-64 pin file
@conda
def test_conda_pin_file():
    run(dpath("test_conda_pin_file"), deployment_method={DeploymentMethod.CONDA})


@conda
def test_conda_python_script():
    run(dpath("test_conda_python_script"), deployment_method={DeploymentMethod.CONDA})


@conda
def test_conda_python_3_7_script():
    run(
        dpath("test_conda_python_3_7_script"),
        deployment_method={DeploymentMethod.CONDA},
    )


@conda
def test_prebuilt_conda_script():
    sp.run(
        f"conda env create -f {dpath('test_prebuilt_conda_script/env.yaml')}",
        shell=True,
    )
    run(dpath("test_prebuilt_conda_script"), deployment_method={DeploymentMethod.CONDA})


@skip_on_windows
@conda
def test_conda_global():
    run(
        dpath("test_conda_global"),
        deployment_method={DeploymentMethod.CONDA},
        executor="dryrun",
    )


@conda
def test_script_pre_py39():
    run(dpath("test_script_pre_py39"), deployment_method={DeploymentMethod.CONDA})


@conda
def test_resource_string_in_cli_or_profile():
    test_path = dpath("test_resource_string_in_cli_or_profile")
    profile = os.path.join(test_path, "profiles")
    # workflow profile is loaded by default
    run(
        test_path,
        snakefile="Snakefile",
        shellcmd=f"snakemake --workflow-profile {profile} -c1 --default-resources slurm_account=foo other_resource='--test'",
    )


@skip_on_windows
@conda
def test_script_xsh():
    run(
        dpath("test_script_xsh"),
        deployment_method={DeploymentMethod.CONDA},
    )


@conda
def test_conda_run():
    run(dpath("test_conda_run"), deployment_method={DeploymentMethod.CONDA})


@conda
def test_issue_3192():
    sp.run(
        "conda create -n test_issue3192 python",
        shell=True,
    )
    run(dpath("test_issue3192"), deployment_method={DeploymentMethod.CONDA})


# Test that container and conda can be run independently using sdm
@skip_on_windows
@apptainer
@connected
@conda
def test_issue_3202():
    run(dpath("test_issue_3202"), deployment_method={DeploymentMethod.APPTAINER})
    run(dpath("test_issue_3202"), deployment_method={DeploymentMethod.CONDA})


# These tests have no explicit dependency on Conda and do not build new conda envs,
# but will fail if 'conda info --json' does not work as expected, because the wrapper
# code uses this to examine the installed software environment.


def test_get_log_none():
    run(dpath("test_get_log_none"))


def test_get_log_both():
    run(dpath("test_get_log_both"))


def test_get_log_stderr():
    run(dpath("test_get_log_stderr"))


def test_get_log_stdout():
    run(dpath("test_get_log_stdout"))


def test_get_log_complex():
    run(dpath("test_get_log_complex"))
