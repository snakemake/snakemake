from snakemake.settings.types import DeploymentSettings

__authors__ = ["Tobias Marschall", "Marcel Martin", "Johannes Köster"]
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import shutil
import signal
import sys
import subprocess as sp
from pathlib import Path
import tempfile
import time

import pytest

sys.path.insert(0, os.path.dirname(__file__))

from .common import run, dpath, apptainer, conda, connected, conda_plugin
from .conftest import (
    skip_on_macos_arm,
    skip_on_windows,
    only_on_windows,
    ON_LINUX,
    ON_WINDOWS,
)

xfail_permissionerror_on_win = (
    pytest.mark.xfail(raises=PermissionError) if ON_WINDOWS else lambda x: x
)


@skip_on_windows
@skip_on_macos_arm
@conda_plugin
def test_script():
    run(
        dpath("test_script"),
        deployment_method={"conda"},
        check_md5=False,
    )


@skip_on_windows
@conda_plugin
def test_script_rs():
    run(
        dpath("test_script_rs"),
        deployment_method={"conda"},
        check_md5=False,
    )


@conda_plugin
def test_conda():
    run(dpath("test_conda"), deployment_method={"conda"})


@conda_plugin
def test_list_software_envs():
    run(dpath("test_conda"), list_software_envs=True, check_results=False)


# TODO failing with FAILED tests/tests.py::test_conda_create_envs_only -
# PermissionError: [WinError 32] The process cannot access the file because
# it is being used by another process:
# 'C:\\Users\\RUNNER~1\\AppData\\Local\\Temp\\snakemake-2q4osog0\\test-env.yaml'
@skip_on_windows
@conda_plugin
def test_cache_or_deploy_software_envs():
    tmpdir = run(
        dpath("test_conda"),
        deployment_method={"conda"},
        cache_or_deploy_software_envs=True,
        check_results=False,
        cleanup=False,
        cleanup_scripts=False,
    )
    assert tmpdir is not None
    env_dir = next(
        (
            p
            for p in (
                Path(tmpdir) / DeploymentSettings().deployment_prefix / "conda"
            ).iterdir()
            if p.is_dir()
        ),
        None,
    )
    assert env_dir is not None
    shutil.rmtree(tmpdir)


@conda_plugin
def test_upstream_conda():
    run(
        dpath("test_conda"),
        deployment_method={"conda"},
    )


@skip_on_windows
@conda_plugin
def test_deploy_script():
    run(dpath("test_deploy_script"), deployment_method={"conda"})


@skip_on_windows
@conda_plugin
def test_deploy_hashing():
    tmpdir = run(
        dpath("test_deploy_hashing"),
        deployment_method={"conda"},
        cleanup=False,
    )
    deployment_prefix = DeploymentSettings().deployment_prefix
    assert tmpdir is not None

    deployed_envs = [
        path
        for path in (Path(tmpdir) / deployment_prefix / "conda").iterdir()
        if path.is_dir() and not path.name == "staged_packages"
    ]
    assert len(deployed_envs) == 1


@conda_plugin
def test_custom_software_deployment_prefix():
    run(
        dpath("test_conda_custom_prefix"),
        deployment_method={"conda"},
        deployment_prefix=".snakemake/custom",
        set_pythonpath=False,
    )


@only_on_windows
@conda_plugin
def test_conda_cmd_exe():
    # Tests the conda environment activation when cmd.exe
    # is used as the shell
    run(dpath("test_conda_cmd_exe"), deployment_method={"conda"})


@skip_on_windows  # wrappers are for linux and macos only
@conda_plugin
def test_wrapper():
    run(
        dpath("test_wrapper"),
        deployment_method={"conda"},
        check_md5=False,
    )


@skip_on_windows  # wrappers are for linux and macos only
@conda_plugin
def test_wrapper_qsub():
    run(
        dpath("test_wrapper"),
        deployment_method={"conda"},
        check_md5=False,
        cluster="./qsub",
    )


@skip_on_windows  # wrappers are for linux and macos only
@conda_plugin
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
            deployment_method={"conda"},
            wrapper_prefix=f"git+file://{tmpdir}",
            check_md5=False,
        )


@skip_on_windows
@connected
@conda_plugin
def test_singularity_conda():
    run(
        dpath("test_singularity_conda"),
        deployment_method={"conda", "container"},
    )


@conda_plugin
@pytest.mark.needs_envmodules
def test_archive():
    run(dpath("test_archive"), archive="workflow-archive.tar.gz")


@skip_on_windows
@conda_plugin
def test_issue635():
    run(
        dpath("test_issue635"),
        deployment_method={"conda"},
        check_md5=False,
    )


@conda_plugin
def test_issue1046():
    run(dpath("test_issue1046"))


@skip_on_windows
@conda_plugin
def test_issue1093():
    run(dpath("test_issue1093"), deployment_method={"conda"})


@conda_plugin
def test_jupyter_notebook():
    run(dpath("test_jupyter_notebook"), deployment_method={"conda"})


@conda_plugin
def test_jupyter_notebook_draft():
    from snakemake.settings.types import NotebookEditMode

    run(
        dpath("test_jupyter_notebook_draft"),
        deployment_method={"conda"},
        edit_notebook=NotebookEditMode(draft_only=True),
        targets=["results/result_intermediate.txt"],
        check_md5=False,
    )


# TODO re-enable as soon as possible:
# The conda env hash has been improved but changed due to that.
# We have to update the containerization test container with the new hashes!
# @skip_on_windows
# @conda_plugin
# def test_containerized():
#     run(
#         dpath("test_containerized"),
#         deployment_method={"conda", "container"},
#     )


@skip_on_windows
@conda_plugin
def test_containerize():
    run(dpath("test_conda"), containerize=True, check_results=False)


@skip_on_windows
@conda_plugin
def test_containerize_wrapper_apptainer():
    tmpdir = None
    try:
        tmpdir = run(
            dpath("test_wrapper"),
            shellcmd="snakemake --containerize apptainer > apptainer.def",
            check_results=False,
            cleanup=False,
            deployment_method={"conda"},
        )
        tmpdir = Path(tmpdir)

        apptainer_def_path = tmpdir / "apptainer.def"
        assert apptainer_def_path.exists(), "apptainer.def was not generated."

        with open(apptainer_def_path) as f:
            apptainer_def_content = f.read()

        assert (
            "#   source: https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v8.1.1/bio/bgzip/environment.yaml"
            in apptainer_def_content
        )

    finally:
        if tmpdir and os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)


def test_converting_path_for_r_script():
    run(
        dpath("test_converting_path_for_r_script"),
        cores=1,
        deployment_method={"conda"},
    )


@skip_on_windows
@conda
def test_conda_named():
    run(dpath("test_conda_named"), deployment_method={"conda"})


# TODO add this again as soon as possible!
# @skip_on_windows
# @conda_plugin
# def test_conda_function():
#     run(
#         dpath("test_conda_function"),
#         deployment_method={"conda"},
#         cores=1,
#     )


@pytest.mark.skipif(not ON_LINUX, reason="This testcase only has a linux-64 pin file")
def test_conda_pin_file():
    run(dpath("test_conda_pin_file"), deployment_method={"conda"})


@conda_plugin
def test_conda_python_script():
    run(dpath("test_conda_python_script"), deployment_method={"conda"})


@skip_on_macos_arm
@conda_plugin
def test_conda_python_3_7_script():
    run(
        dpath("test_conda_python_3_7_script"),
        deployment_method={"conda"},
    )


@conda_plugin
def test_prebuilt_conda_script():
    sp.run(
        f"conda env create -f {dpath('test_prebuilt_conda_script/env.yaml')}",
        shell=True,
    )
    run(dpath("test_prebuilt_conda_script"), deployment_method={"conda"})


@skip_on_windows
@conda_plugin
def test_conda_global():
    run(
        dpath("test_conda_global"),
        deployment_method={"conda"},
        executor="dryrun",
    )


@skip_on_macos_arm
@conda_plugin
def test_script_pre_py39():
    run(dpath("test_script_pre_py39"), deployment_method={"conda"})


@conda_plugin
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
@conda_plugin
def test_script_xsh():
    run(
        dpath("test_script_xsh"),
        deployment_method={"conda"},
    )


@conda_plugin
def test_conda_run():
    run(dpath("test_conda_run"), deployment_method={"conda"})


@conda_plugin
def test_issue_3192():
    with tempfile.TemporaryDirectory() as tmpdir:
        assert (
            sp.run(
                f"conda create -c conda-forge --prefix {tmpdir} python",
                shell=True,
            ).returncode
            == 0
        )
        run(
            dpath("test_issue3192"),
            deployment_method={"conda"},
            config={"env_path": tmpdir},
        )


# Test that container and conda can be run independently using sdm
@skip_on_windows
@connected
@conda_plugin
def test_issue_3202():
    run(dpath("test_issue_3202"), deployment_method={"container"})
    run(dpath("test_issue_3202"), deployment_method={"conda"})


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


@skip_on_windows
@conda_plugin
def test_containerize_checkpoint():
    """Test that containerize considers rules not in the initial DAG (e.g. after checkpoints)."""
    tmpdir = None
    try:
        tmpdir = run(
            dpath("test_containerize_checkpoint"),
            shellcmd="snakemake --containerize > Dockerfile",
            check_results=False,
            cleanup=False,
            deployment_method={"conda"},
        )
        tmpdir = Path(tmpdir)

        dockerfile_path: Path = tmpdir / "Dockerfile"
        assert dockerfile_path.exists(), "Dockerfile was not generated."

        with open(dockerfile_path) as f:
            dockerfile_content = f.read()

        assert "#   source: envs/a.yaml" in dockerfile_content
        assert "#   source: envs/b.yaml" in dockerfile_content
        module_env_path = os.path.join("workflow", "envs", "c.yaml")
        assert f"#   source: {module_env_path}" in dockerfile_content

        # check that COPY/ADD instructions use correct relative paths
        assert "COPY envs/a.yaml" in dockerfile_content
        assert "COPY envs/b.yaml" in dockerfile_content
        assert f"COPY {module_env_path}" in dockerfile_content

        # check three unique environments are being created
        assert (
            dockerfile_content.count("conda env create") == 3
        ), "Expected 3 conda environments to be created."

    finally:
        if tmpdir and os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)


@conda_plugin
def test_issue_1266():
    run(
        dpath("test_github_issue1266"),
        deployment_method={"conda"},
    )


@skip_on_windows
@conda_plugin
def test_github_issue2213():
    run(
        dpath("test_github_issue2213"),
        deployment_method={"conda"},
        check_results=False,
    )


# Parallel deployment tests. Each run has its own working directory, but all
# runs share one software deployment prefix, as on a cluster.
PARALLEL_DEPLOY_TIMEOUT = 900


def _start_parallel_deploy_run(base: Path, name: str, pause=False) -> sp.Popen:
    """Start a run. With pause, the run holds its deployment until the test
    calls _release_deployment (see test_conda_parallel_deploy/run_snakemake.py).
    """
    workdir = base / name
    workdir.mkdir()
    log = open(base / f"{name}.log", "w")
    env = dict(os.environ)
    if pause:
        env["SNAKEMAKE_TEST_DEPLOY_STARTED"] = str(base / "deploy_started")
        env["SNAKEMAKE_TEST_DEPLOY_RELEASE"] = str(base / "deploy_release")
    # A new session, so that the test can kill all child processes.
    return sp.Popen(
        [
            sys.executable,
            str(dpath("test_conda_parallel_deploy/run_snakemake.py")),
            "--snakefile",
            str(dpath("test_conda_parallel_deploy/Snakefile")),
            "--directory",
            str(workdir),
            "--cores",
            "1",
            "--sdm",
            "conda",
            "--sdm-prefix",
            str(base / "envs"),
            "--sdm-cache-prefix",
            str(base / "cache"),
        ],
        stdout=log,
        stderr=sp.STDOUT,
        env=env,
        start_new_session=True,
    )


def _release_deployment(base: Path) -> None:
    (base / "deploy_release").touch()


def _wait_parallel_deploy_run(proc: sp.Popen) -> int:
    try:
        return proc.wait(timeout=PARALLEL_DEPLOY_TIMEOUT)
    except sp.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
        raise


def _parallel_deploy_env_dirs(base: Path):
    prefix = base / "envs" / "conda"
    if not prefix.exists():
        return []
    return [p for p in prefix.iterdir() if p.is_dir() and p.name != "staged_packages"]


def _wait_for_log(base: Path, name: str, proc: sp.Popen, condition, what: str):
    deadline = time.time() + PARALLEL_DEPLOY_TIMEOUT
    while time.time() < deadline:
        if condition():
            return
        if proc.poll() is not None:
            pytest.fail(
                f"Run {name} ended before {what}:\n"
                + (base / f"{name}.log").read_text()
            )
        time.sleep(0.05)
    os.killpg(proc.pid, signal.SIGKILL)
    pytest.fail(f"Run {name} did not reach {what} in time.")


def _wait_for_deployment_start(base: Path, proc: sp.Popen) -> Path:
    """Wait until the paused run A has deployed, but not written the done flag."""
    _wait_for_log(
        base,
        "run_a",
        proc,
        lambda: (base / "deploy_started").exists(),
        "the end of the deployment",
    )
    (env_dir,) = _parallel_deploy_env_dirs(base)
    assert env_dir.with_name(f"{env_dir.name}.lock").exists()
    assert not env_dir.with_name(f"{env_dir.name}.done").exists()
    return env_dir


def _wait_for_lock_wait(base: Path, name: str, proc: sp.Popen) -> None:
    _wait_for_log(
        base,
        name,
        proc,
        lambda: "Waiting for parallel deployment to finish"
        in (base / f"{name}.log").read_text(),
        "the wait for the deployment lock",
    )


def _check_parallel_deploy_run(base: Path, name: str, returncode: int) -> str:
    log = (base / f"{name}.log").read_text()
    assert returncode == 0, log
    assert (base / name / "sum.txt").read_text().strip() == "55"
    assert (base / name / "mean.txt").read_text().strip() == "2.5"
    return log


@skip_on_windows
@connected
@conda_plugin
def test_conda_parallel_deploy(tmp_path):
    run_a = _start_parallel_deploy_run(tmp_path, "run_a", pause=True)
    env_dir = _wait_for_deployment_start(tmp_path, run_a)
    run_b = _start_parallel_deploy_run(tmp_path, "run_b")
    # Run B must wait for the lock of run A, and not use the environment yet.
    _wait_for_lock_wait(tmp_path, "run_b", run_b)
    assert not (tmp_path / "run_b" / "sum.txt").exists()
    _release_deployment(tmp_path)
    rc_a = _wait_parallel_deploy_run(run_a)
    rc_b = _wait_parallel_deploy_run(run_b)
    _check_parallel_deploy_run(tmp_path, "run_a", rc_a)
    _check_parallel_deploy_run(tmp_path, "run_b", rc_b)
    assert env_dir.with_name(f"{env_dir.name}.done").exists()
    assert not env_dir.with_name(f"{env_dir.name}.lock").exists()


@skip_on_windows
@connected
@conda_plugin
def test_conda_killed_deploy(tmp_path):
    # Run A is killed after the deployment, before it writes the done flag.
    run_a = _start_parallel_deploy_run(tmp_path, "run_a", pause=True)
    env_dir = _wait_for_deployment_start(tmp_path, run_a)
    os.killpg(run_a.pid, signal.SIGKILL)
    run_a.wait()
    # The lock of the dead run stays until its lifetime ends.
    assert env_dir.with_name(f"{env_dir.name}.lock").exists()
    assert not env_dir.with_name(f"{env_dir.name}.done").exists()

    run_b = _start_parallel_deploy_run(tmp_path, "run_b")
    log = _check_parallel_deploy_run(
        tmp_path, "run_b", _wait_parallel_deploy_run(run_b)
    )
    assert "Waiting for parallel deployment to finish" in log
    assert "Removing incomplete deployment" in log
    assert env_dir.with_name(f"{env_dir.name}.done").exists()


@skip_on_windows
@connected
@conda_plugin
def test_conda_partial_deploy_is_repaired(tmp_path):
    run_a = _start_parallel_deploy_run(tmp_path, "run_a")
    _check_parallel_deploy_run(tmp_path, "run_a", _wait_parallel_deploy_run(run_a))
    (env_dir,) = _parallel_deploy_env_dirs(tmp_path)

    # Make a partial environment without a lock and without a done flag.
    env_dir.with_name(f"{env_dir.name}.done").unlink()
    (env_dir / "bin" / "datamash").unlink()

    run_b = _start_parallel_deploy_run(tmp_path, "run_b")
    log = _check_parallel_deploy_run(
        tmp_path, "run_b", _wait_parallel_deploy_run(run_b)
    )
    assert "Removing incomplete deployment" in log
    assert (env_dir / "bin" / "datamash").exists()
    assert env_dir.with_name(f"{env_dir.name}.done").exists()
