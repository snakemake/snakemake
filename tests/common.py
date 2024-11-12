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

from snakemake_interface_executor_plugins.settings import SharedFSUsage
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry

from snakemake import api
from snakemake.common import ON_WINDOWS
from snakemake.report.html_reporter import ReportSettings
from snakemake.resources import ResourceScopes
from snakemake.settings import types as settings


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


def has_apptainer():
    return (shutil.which("apptainer") is not None) or (
        shutil.which("singularity") is not None
    )


def has_conda():
    return shutil.which("conda") is not None


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

apptainer = pytest.mark.skipif(
    not has_apptainer(),
    reason="Skipping Apptainer tests because no "
    "apptainer/singularity executable available.",
)

conda = pytest.mark.skipif(
    not has_conda(),
    reason="Skipping Conda tests because no conda executable available.",
)

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
    conda_frontend="conda",
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
    conda_create_envs=False,
    conda_prefix=None,
    wrapper_prefix=None,
    printshellcmds=False,
    default_storage_provider=None,
    default_storage_prefix=None,
    local_storage_prefix=Path(".snakemake/storage"),
    remote_job_local_storage_prefix=None,
    archive=None,
    cluster=None,
    cluster_status=None,
    retries=0,
    resources=dict(),
    default_resources=None,
    group_components=dict(),
    max_threads=None,
    overwrite_groups=dict(),
    configfiles=list(),
    overwrite_resources=dict(),
    batch=None,
    envvars=list(),
    cache=None,
    edit_notebook=None,
    overwrite_scatter=dict(),
    generate_unit_tests=None,
    force_incomplete=False,
    containerize=False,
    forceall=False,
    all_temp=False,
    cleanup_metadata=None,
    rerun_triggers=settings.RerunTrigger.all(),
    storage_provider_settings=None,
    shared_fs_usage=None,
    benchmark_extended=False,
    apptainer_args="",
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
    if check_results:
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
                res = subprocess.run(
                    shellcmd,
                    cwd=path if no_tmpdir else tmpdir,
                    check=True,
                    shell=True,
                    stderr=subprocess.STDOUT,
                    stdout=subprocess.PIPE,
                )
                print(res.stdout.decode())
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
                    if success:
                        print(process.stdout.read().decode())
        except subprocess.CalledProcessError as e:
            success = False
            print(e.stdout.decode(), file=sys.stderr)
    else:
        assert sigint_after is None, "Cannot sent SIGINT when calling directly"

        if cluster is not None:
            executor = "cluster-generic"
            plugin = ExecutorPluginRegistry().get_plugin(executor)
            executor_settings = plugin.settings_cls(
                submit_cmd=cluster, status_cmd=cluster_status
            )
            nodes = 3

        if shared_fs_usage is None:
            shared_fs_usage = SharedFSUsage.all()

        success = True

        with api.SnakemakeApi(
            settings.OutputSettings(
                verbose=True,
                printshellcmds=printshellcmds,
                show_failed_logs=True,
            ),
        ) as snakemake_api:
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
                        overwrite_resources=overwrite_resources,
                        resources=resources,
                        default_resources=default_resources,
                        max_threads=max_threads,
                        overwrite_scatter=overwrite_scatter,
                    ),
                    config_settings=settings.ConfigSettings(
                        config=config,
                        configfiles=configfiles,
                    ),
                    storage_settings=settings.StorageSettings(
                        default_storage_provider=default_storage_provider,
                        default_storage_prefix=default_storage_prefix,
                        all_temp=all_temp,
                        shared_fs_usage=shared_fs_usage,
                        local_storage_prefix=local_storage_prefix,
                        remote_job_local_storage_prefix=remote_job_local_storage_prefix,
                    ),
                    storage_provider_settings=storage_provider_settings,
                    workflow_settings=settings.WorkflowSettings(
                        wrapper_prefix=wrapper_prefix,
                        cache=cache,
                    ),
                    deployment_settings=settings.DeploymentSettings(
                        conda_frontend=conda_frontend,
                        conda_prefix=conda_prefix,
                        deployment_method=deployment_method,
                        apptainer_args=apptainer_args,
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
                        batch=batch,
                        force_incomplete=force_incomplete,
                        forceall=forceall,
                        rerun_triggers=rerun_triggers,
                    ),
                )

                if report is not None:
                    if report_stylesheet is not None:
                        report_stylesheet = Path(report_stylesheet)
                    report_settings = ReportSettings(
                        path=Path(report), stylesheet_path=report_stylesheet
                    )
                    dag_api.create_report(
                        reporter="html",
                        report_settings=report_settings,
                    )
                elif conda_create_envs:
                    dag_api.conda_create_envs()
                elif conda_list_envs:
                    dag_api.conda_list_envs()
                elif archive is not None:
                    dag_api.archive(Path(archive))
                elif generate_unit_tests is not None:
                    dag_api.generate_unit_tests(Path(generate_unit_tests))
                elif containerize:
                    dag_api.containerize()
                elif cleanup_metadata:
                    dag_api.cleanup_metadata(cleanup_metadata)
                else:
                    dag_api.execute_workflow(
                        executor=executor,
                        execution_settings=settings.ExecutionSettings(
                            cleanup_scripts=cleanup_scripts,
                            shadow_prefix=shadow_prefix,
                            retries=retries,
                            edit_notebook=edit_notebook,
                        ),
                        remote_execution_settings=settings.RemoteExecutionSettings(
                            container_image=container_image,
                            seconds_between_status_checks=0,
                            envvars=envvars,
                        ),
                        scheduling_settings=settings.SchedulingSettings(
                            ilp_solver=scheduler_ilp_solver,
                        ),
                        group_settings=settings.GroupSettings(
                            group_components=group_components,
                            overwrite_groups=overwrite_groups,
                        ),
                        executor_settings=executor_settings,
                    )
            except Exception as e:
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

    if not cleanup:
        return tmpdir
    shutil.rmtree(tmpdir, ignore_errors=ON_WINDOWS)
