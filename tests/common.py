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
import urllib.request
import pytest
import glob
import subprocess
import tarfile
from typing import TypeAlias

from snakemake_interface_executor_plugins.settings import SharedFSUsage
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry

from snakemake import api
from snakemake.common import ON_WINDOWS
from snakemake.report.html_reporter import ReportSettings
from snakemake.resources import ResourceScopes
from snakemake.scheduling.milp import SchedulerSettings
from snakemake.settings import types as settings


#: File system path as string or pathlike object.
StrPath: TypeAlias = str | os.PathLike


def dpath(path: StrPath) -> Path:
    """get the path to a data file (relative to the directory this
    test lives in)"""
    return (Path(__file__).parent / path).resolve()


def md5sum(filename: StrPath, ignore_newlines: bool = False) -> str:
    if ignore_newlines:
        with open(filename, "r", encoding="utf-8", errors="surrogateescape") as f:
            data = f.read().strip().encode("utf8", errors="surrogateescape")
    else:
        data = open(filename, "rb").read().strip()
    return hashlib.md5(data, usedforsecurity=False).hexdigest()


# test skipping
def is_connected() -> bool:
    if is_ci():
        # always assume internet connection in CI
        return True
    try:
        urllib.request.urlopen("https://www.google.com", timeout=1)
        return True
    except (urllib.request.URLError, TimeoutError):
        return False


def is_ci() -> bool:
    return "CI" in os.environ


def has_gcloud_service_key() -> bool:
    return "GCP_AVAILABLE" in os.environ


def has_azbatch_account_url() -> str | None:
    return os.environ.get("AZ_BATCH_ACCOUNT_URL")


def has_zenodo_token() -> str | None:
    return os.environ.get("ZENODO_SANDBOX_PAT")


def has_apptainer() -> bool:
    return (shutil.which("apptainer") is not None) or (
        shutil.which("singularity") is not None
    )


def has_conda() -> bool:
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


def copy(src: StrPath, dst: StrPath):
    if os.path.isdir(src):
        shutil.copytree(src, os.path.join(dst, os.path.basename(src)))
    else:
        shutil.copy(src, dst)


def get_expected_files(results_dir: StrPath) -> list[str]:
    """Recursively walk through the expected-results directory to enumerate
    all expected files."""
    return [
        os.path.relpath(f, results_dir)
        for f in glob.iglob(os.path.join(results_dir, "**/**"), recursive=True)
        if not os.path.isdir(f)
    ]


def untar_folder(tar_file: StrPath, output_path: StrPath):
    if not os.path.isdir(output_path):
        with tarfile.open(tar_file) as tar:
            tar.extractall(path=output_path)


def print_tree(path: str, exclude: StrPath | None = None):
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
    path: StrPath,
    shouldfail: bool = False,
    snakefile: StrPath = "Snakefile",
    no_tmpdir: bool = False,
    check_md5: bool = True,
    check_results: bool | None = None,
    cores: int = 3,
    nodes: int | None = None,
    set_pythonpath: bool = True,
    cleanup: bool = True,
    conda_frontend="conda",
    config=dict(),
    targets=set(),
    container_image=os.environ.get("CONTAINER_IMAGE", "snakemake/snakemake:latest"),
    shellcmd: str | None = None,
    sigint_after: float | None = None,
    overwrite_resource_scopes=None,
    executor="local",
    executor_settings=None,
    cleanup_scripts=True,
    scheduler_ilp_solver=None,
    report=None,
    report_after_run=False,
    report_stylesheet=None,
    report_metadata=None,
    deployment_method=frozenset(),
    shadow_prefix=None,
    until=frozenset(),
    omit_from=frozenset(),
    forcerun=frozenset(),
    trust_io_cache=False,
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
    tmpdir: StrPath | None = None,
) -> Path | None:
    """
    Test the Snakefile in the path.
    There must be a Snakefile in the path and a subdirectory named
    expected-results. If cleanup is False, we return the temporary
    directory to the calling test for inspection, and the test should
    clean it up.

    Parameters
    ----------
    path
        Path containing workflow to run.
    shouldfail
        Whether the run is expected to fail.
    snakefile
        Path to Snakefile, relative to ``path``.
    shellcmd
        Shell command to run. Must start with "snakemake". If given, Snakemake will be run in a
        subprocess.
    sigint_after
        If not None, send a SIGINT signal after this many seconds.
    tmpdir
        Temporary directory to run in. If None one will be created automatically.
    no_tmpdir
        If true run directly in ``path`` instead of a temporary directory.
    cleanup
        Whether to delete the temporary directory after running.
    set_pythonpath
        If true set the ``PYTHONPATH`` environment variable to current working directory. Otherwise
        ensure it is not set.

    Returns
    -------
    Path | None
        Path to temporary directory if ``cleanup`` is false, otherwise None.
    """
    path = Path(path)

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

    results_dir = path / "expected-results"
    original_snakefile = path / snakefile
    original_dirname = original_snakefile.parent.name
    assert original_snakefile.exists()
    if check_results:
        assert (
            results_dir.exists() and results_dir.is_dir()
        ), f"{results_dir} does not exist"

    if tmpdir is None:
        # If we need to further check results, we won't cleanup tmpdir
        tmpdir = next(tempfile._get_candidate_names())
        tmpdir = os.path.join(
            tempfile.gettempdir(), f"snakemake-{original_dirname}-{tmpdir}"
        )
        os.mkdir(tmpdir)

        # copy files
        for f in os.listdir(path):
            copy(os.path.join(path, f), tmpdir)

    else:
        tmpdir = os.fsdecode(tmpdir)

    # Snakefile is now in temporary directory
    snakefile = join(tmpdir, snakefile)

    snakemake_api = None
    exception = None

    config = dict(config)

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
                        trust_io_cache=trust_io_cache,
                    ),
                )

                if report is not None and not report_after_run:
                    if report_stylesheet is not None:
                        report_stylesheet = Path(report_stylesheet)
                    if report_metadata is not None:
                        report_metadata = Path(report_metadata)
                    report_settings = ReportSettings(
                        path=Path(report), stylesheet_path=report_stylesheet
                    )
                    global_report_settings = settings.GlobalReportSettings(
                        metadata_template=report_metadata
                    )
                    dag_api.create_report(
                        reporter="html",
                        report_settings=report_settings,
                        global_report_settings=global_report_settings,
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
                        group_settings=settings.GroupSettings(
                            group_components=group_components,
                            overwrite_groups=overwrite_groups,
                        ),
                        executor_settings=executor_settings,
                        scheduler_settings=SchedulerSettings(
                            solver=scheduler_ilp_solver,
                        ),
                    )

                if report_after_run and report:
                    if report_stylesheet is not None:
                        report_stylesheet = Path(report_stylesheet)
                    report_settings = ReportSettings(
                        path=Path(report), stylesheet_path=report_stylesheet
                    )
                    dag_api.create_report(
                        reporter="html",
                        report_settings=report_settings,
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
            if exception is not None:
                raise exception
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
        return Path(tmpdir)

    shutil.rmtree(tmpdir, ignore_errors=ON_WINDOWS)
    return None
