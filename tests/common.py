__authors__ = ["Tobias Marschall", "Marcel Martin", "Johannes Köster"]
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
from pathlib import Path
import signal
import sys
import shlex
import shutil
import time
from os.path import join
import tempfile
import hashlib
import urllib
import pytest
import glob
import subprocess
import tarfile
from snakemake import api, settings

from snakemake.common import ON_WINDOWS
from snakemake.resources import ResourceScopes


def dpath(path):
    """get the path to a data file (relative to the directory this
    test lives in)"""
    return os.path.realpath(join(os.path.dirname(__file__), path))


def md5sum(filename, ignore_newlines=False):
    if ignore_newlines:
        with open(filename, "r", encoding="utf-8", errors="surrogateescape") as f:
            data = f.read().strip().encode("utf8", errors="surrogateescape")
    else:
        data = open(filename, "rb").read().strip()
    return hashlib.md5(data).hexdigest()


# test skipping
def is_connected():
    try:
        urllib.request.urlopen("http://www.google.com", timeout=1)
        return True
    except urllib.request.URLError:
        return False


def is_ci():
    return "CI" in os.environ


def has_gcloud_service_key():
    return "GCP_AVAILABLE" in os.environ


def has_azbatch_account_url():
    return os.environ.get("AZ_BATCH_ACCOUNT_URL")


def has_zenodo_token():
    return os.environ.get("ZENODO_SANDBOX_PAT")


gcloud = pytest.mark.skipif(
    not is_connected() or not has_gcloud_service_key(),
    reason="Skipping GCLOUD tests because not on "
    "CI, no inet connection or not logged "
    "in to gcloud.",
)

azbatch = pytest.mark.skipif(
    not is_connected() or not has_azbatch_account_url(),
    reason="Skipping AZBATCH tests because "
    "no inet connection or no AZ_BATCH_ACCOUNT_URL.",
)

connected = pytest.mark.skipif(not is_connected(), reason="no internet connection")

ci = pytest.mark.skipif(not is_ci(), reason="not in CI")
not_ci = pytest.mark.skipif(is_ci(), reason="skipped in CI")

zenodo = pytest.mark.skipif(
    not has_zenodo_token(), reason="no ZENODO_SANDBOX_PAT provided"
)


def copy(src, dst):
    if os.path.isdir(src):
        shutil.copytree(src, os.path.join(dst, os.path.basename(src)))
    else:
        shutil.copy(src, dst)


def get_expected_files(results_dir):
    """Recursively walk through the expected-results directory to enumerate
    all expected files."""
    return [
        os.path.relpath(f, results_dir)
        for f in glob.iglob(os.path.join(results_dir, "**/**"), recursive=True)
        if not os.path.isdir(f)
    ]


def untar_folder(tar_file, output_path):
    if not os.path.isdir(output_path):
        with tarfile.open(tar_file) as tar:
            tar.extractall(path=output_path)


def print_tree(path, exclude=None):
    for root, _dirs, files in os.walk(path):
        if exclude and root.startswith(os.path.join(path, exclude)):
            continue
        level = root.replace(path, "").count(os.sep)
        indent = " " * 4 * level
        print(f"{indent}{os.path.basename(root)}/")
        subindent = " " * 4 * (level + 1)
        for f in files:
            print(f"{subindent}{f}")


def run(
    path,
    shouldfail=False,
    snakefile="Snakefile",
    subpath=None,
    no_tmpdir=False,
    check_md5=True,
    check_results=None,
    cores=3,
    nodes=None,
    set_pythonpath=True,
    cleanup=True,
    conda_frontend="mamba",
    config=dict(),
    targets=set(),
    container_image=os.environ.get("CONTAINER_IMAGE", "snakemake/snakemake:latest"),
    shellcmd=None,
    sigint_after=None,
    overwrite_resource_scopes=None,
    executor="local",
    executor_settings=None,
    cleanup_scripts=True,
    scheduler_ilp_solver=None,
    report=None,
    report_stylesheet=None,
    deployment_method=frozenset(),
    shadow_prefix=None,
    until=frozenset(),
    omit_from=frozenset(),
    forcerun=frozenset(),
    conda_list_envs=False,
    conda_prefix=None,
    wrapper_prefix=None,
    printshellcmds=False,
):
    """
    Test the Snakefile in the path.
    There must be a Snakefile in the path and a subdirectory named
    expected-results. If cleanup is False, we return the temporary
    directory to the calling test for inspection, and the test should
    clean it up.
    """
    if check_results is None:
        if not shouldfail:
            check_results = True
        else:
            check_results = False

    if set_pythonpath:
        # Enforce current workdir (the snakemake source dir) to also be in PYTHONPATH
        # when subprocesses are invoked in the tempdir defined below.
        os.environ["PYTHONPATH"] = os.getcwd()
    elif "PYTHONPATH" in os.environ:
        del os.environ["PYTHONPATH"]

    results_dir = join(path, "expected-results")
    original_snakefile = join(path, snakefile)
    assert os.path.exists(original_snakefile)
    assert os.path.exists(results_dir) and os.path.isdir(
        results_dir
    ), "{} does not exist".format(results_dir)

    # If we need to further check results, we won't cleanup tmpdir
    tmpdir = next(tempfile._get_candidate_names())
    tmpdir = os.path.join(tempfile.gettempdir(), "snakemake-%s" % tmpdir)
    os.mkdir(tmpdir)

    config = dict(config)

    # copy files
    for f in os.listdir(path):
        copy(os.path.join(path, f), tmpdir)

    # Snakefile is now in temporary directory
    snakefile = join(tmpdir, snakefile)

    snakemake_api = None
    exception = None

    # run snakemake
    if shellcmd:
        if not shellcmd.startswith("snakemake"):
            raise ValueError("shellcmd does not start with snakemake")
        shellcmd = "{} -m {}".format(sys.executable, shellcmd)
        try:
            if sigint_after is None:
                subprocess.run(
                    shellcmd,
                    cwd=path if no_tmpdir else tmpdir,
                    check=True,
                    shell=True,
                    stderr=subprocess.STDOUT,
                    stdout=subprocess.PIPE,
                )
                success = True
            else:
                with subprocess.Popen(
                    shlex.split(shellcmd),
                    cwd=path if no_tmpdir else tmpdir,
                    stderr=subprocess.STDOUT,
                    stdout=subprocess.PIPE,
                ) as process:
                    time.sleep(sigint_after)
                    process.send_signal(signal.SIGINT)
                    time.sleep(2)
                    success = process.returncode == 0
        except subprocess.CalledProcessError as e:
            success = False
            print(e.stdout.decode(), file=sys.stderr)
    else:
        assert sigint_after is None, "Cannot sent SIGINT when calling directly"

        success = True

        snakemake_api = api.SnakemakeApi(
            settings.OutputSettings(
                verbose=True,
                printshellcmds=printshellcmds,
            ),
        )

        try:
            workflow_api = snakemake_api.workflow(
                resource_settings=settings.ResourceSettings(
                    cores=cores,
                    nodes=nodes,
                    overwrite_resource_scopes=(
                        ResourceScopes(overwrite_resource_scopes)
                        if overwrite_resource_scopes is not None
                        else dict()
                    ),
                ),
                config_settings=settings.ConfigSettings(
                    config=config,
                ),
                storage_settings=settings.StorageSettings(),
                workflow_settings=settings.WorkflowSettings(
                    wrapper_prefix=wrapper_prefix,
                ),
                snakefile=Path(original_snakefile if no_tmpdir else snakefile),
                workdir=Path(path if no_tmpdir else tmpdir),
            )

            dag_api = workflow_api.dag(
                dag_settings=settings.DAGSettings(
                    targets=targets,
                    until=until,
                    omit_from=omit_from,
                    forcerun=forcerun,
                ),
                deployment_settings=settings.DeploymentSettings(
                    conda_frontend=conda_frontend,
                    conda_prefix=conda_prefix,
                    deployment_method=deployment_method,
                ),
            )

            if report is not None:
                dag_api.create_report(path=report, stylesheet=report_stylesheet)
            elif conda_list_envs:
                dag_api.conda_list_envs()
            else:
                dag_api.execute_workflow(
                    executor=executor,
                    execution_settings=settings.ExecutionSettings(
                        cleanup_scripts=cleanup_scripts,
                        shadow_prefix=shadow_prefix,
                    ),
                    remote_execution_settings=settings.RemoteExecutionSettings(
                        container_image=container_image,
                        seconds_between_status_checks=0,
                    ),
                    scheduling_settings=settings.SchedulingSettings(
                        ilp_solver=scheduler_ilp_solver,
                    ),
                    executor_settings=executor_settings,
                )
        except Exception as e:
            import pdb; pdb.set_trace()
            success = False
            exception = e

    if shouldfail:
        assert not success, "expected error on execution"
    else:
        if not success:
            if snakemake_api is not None and exception is not None:
                snakemake_api.print_exception(exception)
            print("Workdir:")
            print_tree(tmpdir, exclude=".snakemake/conda")
        assert success, "expected successful execution"

    if check_results:
        for resultfile in get_expected_files(results_dir):
            if resultfile in [".gitignore", ".gitkeep"] or not os.path.isfile(
                os.path.join(results_dir, resultfile)
            ):
                # this means tests cannot use directories as output files
                continue
            targetfile = join(tmpdir, resultfile)
            expectedfile = join(results_dir, resultfile)

            if ON_WINDOWS:
                if os.path.exists(join(results_dir, resultfile + "_WIN")):
                    continue  # Skip test if a Windows specific file exists
                if resultfile.endswith("_WIN"):
                    targetfile = join(tmpdir, resultfile[:-4])
            elif resultfile.endswith("_WIN"):
                # Skip win specific result files on Posix platforms
                continue

            assert os.path.exists(targetfile), 'expected file "{}" not produced'.format(
                resultfile
            )
            if check_md5:
                md5expected = md5sum(expectedfile, ignore_newlines=ON_WINDOWS)
                md5target = md5sum(targetfile, ignore_newlines=ON_WINDOWS)
                if md5target != md5expected:
                    with open(expectedfile) as expected:
                        expected_content = expected.read().strip()
                    with open(targetfile) as target:
                        content = target.read().strip()
                    assert (
                        False
                    ), "wrong result produced for file '{resultfile}':\n------found------\n{content}\n-----expected-----\n{expected_content}\n-----------------".format(
                        resultfile=resultfile,
                        content=content,
                        expected_content=expected_content,
                    )

    if snakemake_api is not None:
        snakemake_api.cleanup()

    if not cleanup:
        return tmpdir
    shutil.rmtree(tmpdir, ignore_errors=ON_WINDOWS)
