__authors__ = ["Tobias Marschall", "Marcel Martin", "Johannes Köster"]
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import sys
import uuid
import subprocess as sp
from pathlib import Path

from snakemake import parse_cores_jobs
from snakemake.exceptions import CliException
from snakemake.utils import available_cpu_count

sys.path.insert(0, os.path.dirname(__file__))

from .common import *
from .conftest import skip_on_windows, only_on_windows, ON_WINDOWS, needs_strace


def test_list_untracked():
    run(dpath("test_list_untracked"))


xfail_permissionerror_on_win = (
    pytest.mark.xfail(raises=PermissionError) if ON_WINDOWS else lambda x: x
)


def test_delete_all_output():
    run(dpath("test_delete_all_output"))


def test_github_issue_14():
    """Add cleanup_scripts argument to allow the user to keep scripts"""
    # Return temporary directory for inspection - we should keep scripts here
    tmpdir = run(dpath("test_github_issue_14"), cleanup=False, cleanup_scripts=False)
    assert os.listdir(os.path.join(tmpdir, ".snakemake", "scripts"))
    shutil.rmtree(tmpdir)

    # And not here
    tmpdir = run(dpath("test_github_issue_14"), cleanup=False)
    assert not os.listdir(os.path.join(tmpdir, ".snakemake", "scripts"))
    shutil.rmtree(tmpdir)


def test_issue956():
    run(dpath("test_issue956"))


@skip_on_windows
def test01():
    run(dpath("test01"))


def test02():
    run(dpath("test02"))


def test03():
    run(dpath("test03"), targets=["test.out"])


def test04():
    run(dpath("test04"), targets=["test.out"])


def test05():
    run(dpath("test05"))


def test06():
    run(dpath("test06"), targets=["test.bla.out"])


def test07():
    run(dpath("test07"), targets=["test.out", "test2.out"])


def test08():
    run(dpath("test08"), targets=["test.out", "test2.out"])


def test09():
    run(dpath("test09"), shouldfail=True)


def test10():
    run(dpath("test10"))


def test11():
    run(dpath("test11"))


def test12():
    run(dpath("test12"))


def test13():
    run(dpath("test13"))


@skip_on_windows
def test14():
    os.environ["TESTVAR"] = "test"
    os.environ["TESTVAR2"] = "test"
    run(dpath("test14"), snakefile="Snakefile.nonstandard", cluster="./qsub")


@skip_on_windows
def test_cluster_statusscript():
    os.environ["TESTVAR"] = "test"
    run(
        dpath("test_cluster_statusscript"),
        snakefile="Snakefile.nonstandard",
        cluster="./qsub",
        cluster_status="./status.sh",
    )


@skip_on_windows
def test_cluster_cancelscript():
    outdir = run(
        dpath("test_cluster_cancelscript"),
        snakefile="Snakefile.nonstandard",
        shellcmd=(
            "snakemake -j 10 --cluster=./sbatch --cluster-cancel=./scancel.sh "
            "--cluster-status=./status.sh -s Snakefile.nonstandard"
        ),
        shouldfail=True,
        cleanup=False,
        sigint_after=4,
    )
    scancel_txt = open("%s/scancel.txt" % outdir).read()
    scancel_lines = scancel_txt.splitlines()
    assert len(scancel_lines) == 1
    assert scancel_lines[0].startswith("cancel")
    assert len(scancel_lines[0].split(" ")) == 3


@skip_on_windows
def test_cluster_sidecar():
    run(
        dpath("test_cluster_sidecar"),
        shellcmd=("snakemake -j 10 --cluster=./sbatch --cluster-sidecar=./sidecar.sh"),
    )


@skip_on_windows
def test_cluster_cancelscript_nargs1():
    outdir = run(
        dpath("test_cluster_cancelscript"),
        snakefile="Snakefile.nonstandard",
        shellcmd=(
            "snakemake -j 10 --cluster=./sbatch --cluster-cancel=./scancel.sh "
            "--cluster-status=./status.sh --cluster-cancel-nargs=1 "
            "-s Snakefile.nonstandard"
        ),
        shouldfail=True,
        cleanup=False,
        sigint_after=4,
    )
    scancel_txt = open("%s/scancel.txt" % outdir).read()
    scancel_lines = scancel_txt.splitlines()
    assert len(scancel_lines) == 2
    assert scancel_lines[0].startswith("cancel")
    assert scancel_lines[1].startswith("cancel")
    assert len(scancel_lines[0].split(" ")) == 2
    assert len(scancel_lines[1].split(" ")) == 2


@skip_on_windows
def test_cluster_statusscript_multi():
    os.environ["TESTVAR"] = "test"
    run(
        dpath("test_cluster_statusscript_multi"),
        snakefile="Snakefile.nonstandard",
        cluster="./sbatch",
        cluster_status="./status.sh",
    )


def test15():
    run(dpath("test15"))


def test_glpk_solver():
    run(dpath("test_solver"), scheduler_ilp_solver="GLPK_CMD")


@skip_on_windows
def test_coin_solver():
    run(dpath("test_solver"), scheduler_ilp_solver="COIN_CMD")


def test_directory():
    run(
        dpath("test_directory"),
        targets=[
            "downstream",
            "symlinked_input",
            "child_to_input",
            "not_child_to_other",
        ],
    )
    run(dpath("test_directory"), targets=["file_expecting_dir"], shouldfail=True)
    run(dpath("test_directory"), targets=["dir_expecting_file"], shouldfail=True)
    run(dpath("test_directory"), targets=["child_to_other"], shouldfail=True)


@skip_on_windows
def test_directory2():
    run(
        dpath("test_directory"),
        targets=[
            "downstream",
            "symlinked_input",
            "child_to_input",
            "some/dir-child",
            "some/shadow",
        ],
    )


def test_ancient():
    run(dpath("test_ancient"), targets=["D", "C", "old_file"])


def test_report():
    run(
        dpath("test_report"),
        report="report.html",
        report_stylesheet="custom-stylesheet.css",
        check_md5=False,
    )


def test_report_zip():
    run(dpath("test_report_zip"), report="report.zip", check_md5=False)


def test_report_dir():
    run(dpath("test_report_dir"), report="report.zip", check_md5=False)


def test_report_display_code():
    run(dpath("test_report_display_code"), report="report.html", check_md5=False)


def test_dynamic():
    run(dpath("test_dynamic"))


def test_params():
    run(dpath("test_params"))


def test_same_wildcard():
    run(dpath("test_same_wildcard"))


def test_conditional():
    run(
        dpath("test_conditional"),
        targets="test.out test.0.out test.1.out test.2.out".split(),
    )


def test_unpack_dict():
    run(dpath("test_unpack_dict"))


def test_unpack_list():
    run(dpath("test_unpack_list"))


def test_shell():
    run(dpath("test_shell"))


@skip_on_windows
def test_temp():
    run(dpath("test_temp"), cluster="./qsub", targets="test.realigned.bam".split())


def test_keyword_list():
    run(dpath("test_keyword_list"))


# Fails on WIN because some snakemake doesn't release the logfile
# which cause a PermissionError when the test setup tries to
# remove the temporary files
@skip_on_windows
def test_subworkflows():
    run(dpath("test_subworkflows"), subpath=dpath("test02"))


def test_globwildcards():
    run(dpath("test_globwildcards"))


def test_local_import():
    run(dpath("test_local_import"))


def test_ruledeps():
    run(dpath("test_ruledeps"))


def test_persistent_dict():
    try:
        import pytools

        run(dpath("test_persistent_dict"))
    except ImportError:
        pass


@connected
def test_url_include():
    run(dpath("test_url_include"))


def test_touch():
    run(dpath("test_touch"))


def test_touch_flag_with_directories():
    run(dpath("test_touch_with_directories"), touch=True)


def test_config():
    run(dpath("test_config"))


def test_update_config():
    run(dpath("test_update_config"))


def test_config_merging():
    run(
        dpath("test_config_merging"),
        shellcmd='snakemake -j 1 --configfile config_cmdline_01.yaml config_cmdline_02.yaml --config "block={bowser: cmdline_bowser}" "block={toad: cmdline_toad}"',
    )


def test_wildcard_keyword():
    run(dpath("test_wildcard_keyword"))


@skip_on_windows
def test_benchmark():
    run(dpath("test_benchmark"), check_md5=False)


def test_temp_expand():
    run(dpath("test_temp_expand"))


def test_wildcard_count_ambiguity():
    run(dpath("test_wildcard_count_ambiguity"))


@skip_on_windows
def test_srcdir():
    run(dpath("test_srcdir"))


def test_multiple_includes():
    run(dpath("test_multiple_includes"))


def test_name_override():
    run(dpath("test_name_override"))


def test_yaml_config():
    run(dpath("test_yaml_config"))


@skip_on_windows
@pytest.mark.xfail(
    reason="moto currently fails with \"'_patch' object has no attribute 'is_local'\""
)
def test_remote():
    run(dpath("test_remote"), cores=1)


@skip_on_windows
def test_cluster_sync():
    os.environ["TESTVAR"] = "test"
    run(dpath("test14"), snakefile="Snakefile.nonstandard", cluster_sync="./qsub")


@pytest.mark.skip(reason="This does not work reliably in CircleCI.")
def test_symlink_temp():
    run(dpath("test_symlink_temp"), shouldfail=True)


def test_empty_include():
    run(dpath("test_empty_include"))


@skip_on_windows
def test_script():
    run(dpath("test_script"), use_conda=True, check_md5=False)


def test_script_python():
    run(dpath("test_script_py"))


@skip_on_windows  # Test relies on perl
def test_shadow():
    run(dpath("test_shadow"))


@skip_on_windows
def test_shadow_copy():
    run(dpath("test_shadow_copy"))


@skip_on_windows  # Symbolic link privileges needed to work
def test_shadow_prefix():
    run(dpath("test_shadow_prefix"), shadow_prefix="shadowdir")


@skip_on_windows
def test_shadow_prefix_qsub():
    run(dpath("test_shadow_prefix"), shadow_prefix="shadowdir", cluster="./qsub")


@skip_on_windows
def test_shadowed_log():
    run(dpath("test_shadowed_log"))


def test_until():
    run(
        dpath("test_until"),
        until=[
            "leveltwo_first",  # rule name
            "leveltwo_second.txt",  # file name
            "second_wildcard",
        ],
    )  # wildcard rule


def test_omitfrom():
    run(
        dpath("test_omitfrom"),
        omit_from=[
            "leveltwo_first",  # rule name
            "leveltwo_second.txt",  # file name
            "second_wildcard",
        ],
    )  # wildcard rule


def test_nonstr_params():
    run(dpath("test_nonstr_params"))


def test_delete_output():
    run(dpath("test_delete_output"), cores=1)


def test_input_generator():
    run(dpath("test_input_generator"))


def test_symlink_time_handling():
    # See Snakefile for notes on why this fails on some systems
    if os.utime in os.supports_follow_symlinks:
        run(dpath("test_symlink_time_handling"))


@skip_on_windows
def test_protected_symlink_output():
    run(dpath("test_protected_symlink_output"))


def test_issue328():
    try:
        import pytools

        run(dpath("test_issue328"), forcerun=["split"])
    except ImportError:
        # skip test if import fails
        pass


def test_conda():
    run(dpath("test_conda"), use_conda=True)


def test_conda_list_envs():
    run(dpath("test_conda"), list_conda_envs=True, check_results=False)


def test_upstream_conda():
    run(dpath("test_conda"), use_conda=True, conda_frontend="conda")


@skip_on_windows
def test_deploy_script():
    run(dpath("test_deploy_script"), use_conda=True)


@skip_on_windows
def test_deploy_hashing():
    tmpdir = run(dpath("test_deploy_hashing"), use_conda=True, cleanup=False)
    assert len(next(os.walk(os.path.join(tmpdir, ".snakemake/conda")))[1]) == 2


def test_conda_custom_prefix():
    run(
        dpath("test_conda_custom_prefix"),
        use_conda=True,
        conda_prefix="custom",
        set_pythonpath=False,
    )


@only_on_windows
def test_conda_cmd_exe():
    # Tests the conda environment activation when cmd.exe
    # is used as the shell
    run(dpath("test_conda_cmd_exe"), use_conda=True)


@skip_on_windows  # wrappers are for linux and macos only
def test_wrapper():
    run(dpath("test_wrapper"), use_conda=True)


@skip_on_windows  # wrappers are for linux and macos only
def test_wrapper_local_git_prefix():
    import git

    with tempfile.TemporaryDirectory() as tmpdir:
        print("Cloning wrapper repo...")
        repo = git.Repo.clone_from(
            "https://github.com/snakemake/snakemake-wrappers", tmpdir
        )
        print("Cloning complete.")

        run(
            dpath("test_wrapper"), use_conda=True, wrapper_prefix=f"git+file://{tmpdir}"
        )


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


def test_spaces_in_fnames():
    run(
        dpath("test_spaces_in_fnames"),
        # cluster="./qsub",
        targets=["test bam file realigned.bam"],
        printshellcmds=True,
    )


# TODO deactivate because of problems with moto and boto3.
# def test_static_remote():
#     import importlib
#     try:
#         importlib.reload(boto3)
#         importlib.reload(moto)
#         # only run the remote file test if the dependencies
#         # are installed, otherwise do nothing
#         run(dpath("test_static_remote"), cores=1)
#     except ImportError:
#         pass


@connected
def test_remote_ncbi_simple():
    try:
        import Bio

        # only run the remote file test if the dependencies
        # are installed, otherwise do nothing
        run(dpath("test_remote_ncbi_simple"))
    except ImportError:
        pass


@connected
def test_remote_ncbi():
    try:
        import Bio

        # only run the remote file test if the dependencies
        # are installed, otherwise do nothing
        run(dpath("test_remote_ncbi"))
    except ImportError:
        pass


@ci
@skip_on_windows
def test_remote_irods():
    run(dpath("test_remote_irods"))


def test_deferred_func_eval():
    run(dpath("test_deferred_func_eval"))


def test_format_params():
    run(dpath("test_format_params"), check_md5=True)


def test_rule_defined_in_for_loop():
    # issue 257
    run(dpath("test_rule_defined_in_for_loop"))


def test_issue381():
    run(dpath("test_issue381"))


def test_format_wildcards():
    run(dpath("test_format_wildcards"))


def test_with_parentheses():
    run(dpath("test (with parenthese's)"))


def test_dup_out_patterns():
    """Duplicate output patterns should emit an error

    Duplicate output patterns can be detected on the rule level
    """
    run(dpath("test_dup_out_patterns"), shouldfail=True)


@skip_on_windows
def test_restartable_job_cmd_exit_1_no_restart():
    """Test the restartable job feature on ``exit 1``

    The shell snippet in the Snakemake file will fail the first time
    and succeed the second time.
    """
    run(
        dpath("test_restartable_job_cmd_exit_1"),
        cluster="./qsub",
        restart_times=0,
        shouldfail=True,
    )


@skip_on_windows
def test_restartable_job_cmd_exit_1_one_restart():
    # Restarting once is enough
    run(
        dpath("test_restartable_job_cmd_exit_1"),
        cluster="./qsub",
        restart_times=1,
        printshellcmds=True,
    )


@skip_on_windows
def test_restartable_job_qsub_exit_1():
    """Test the restartable job feature when qsub fails

    The qsub in the subdirectory will fail the first time and succeed the
    second time.
    """
    # Even two consecutive times should fail as files are cleared
    run(
        dpath("test_restartable_job_qsub_exit_1"),
        cluster="./qsub",
        restart_times=0,
        shouldfail=True,
    )
    run(
        dpath("test_restartable_job_qsub_exit_1"),
        cluster="./qsub",
        restart_times=0,
        shouldfail=True,
    )
    # Restarting once is enough
    run(
        dpath("test_restartable_job_qsub_exit_1"),
        cluster="./qsub",
        restart_times=1,
        shouldfail=False,
    )


def test_threads():
    run(dpath("test_threads"), cores=20)


def test_threads0():
    run(dpath("test_threads0"))


def test_dynamic_temp():
    run(dpath("test_dynamic_temp"))


# TODO this currently hangs. Has to be investigated (issue #660).
# def test_ftp_immediate_close():
#    try:
#        import ftputil
#
#        # only run the remote file test if the dependencies
#        # are installed, otherwise do nothing
#        run(dpath("test_ftp_immediate_close"))
#    except ImportError:
#        pass


def test_issue260():
    run(dpath("test_issue260"))


@skip_on_windows
@not_ci
def test_default_remote():
    run(
        dpath("test_default_remote"),
        cores=1,
        default_remote_provider="S3Mocked",
        default_remote_prefix="test-remote-bucket",
    )


def test_run_namedlist():
    run(dpath("test_run_namedlist"))


@connected
@not_ci
@skip_on_windows
def test_remote_gs():
    run(dpath("test_remote_gs"))


@pytest.mark.skip(reason="Need to choose how to provide billable project")
@connected
@not_ci
def test_gs_requester_pays(
    requesting_project=None,
    requesting_url="gcp-public-data-landsat/LC08/01/001/003/LC08_L1GT_001003_20170430_20170501_01_RT/LC08_L1GT_001003_20170430_20170501_01_RT_MTL.txt",
):
    """Tests pull-request 79 / issue 96 for billable user projects on GS

    If requesting_project None, behaves as test_remote_gs().

    Parameters
    ----------
    requesting_project: Optional[str]
        User project to bill for download. None will not provide the project for
        requester-pays as is the usual default
    requesting_url: str
        URL of the bucket to download. The default will match the expected output but is a
        bucket that doesn't require requester pays.
    """
    # create temporary config file
    with tempfile.NamedTemporaryFile(suffix=".yaml") as handle:
        # specify project and url for download
        if requesting_project is None:
            handle.write(b"project: null\n")
        else:
            handle.write('project: "{}"\n'.format(requesting_project).encode())
        handle.write('url: "{}"\n'.format(requesting_url).encode())
        # make sure we can read them
        handle.flush()
        # run the pipeline
        run(dpath("test_gs_requester_pays"), configfiles=[handle.name], forceall=True)


@pytest.mark.skip(reason="We need free azure access to test this in CircleCI.")
@connected
@ci
@skip_on_windows
def test_remote_azure():
    run(dpath("test_remote_azure"))


def test_remote_log():
    run(dpath("test_remote_log"), shouldfail=True)


@connected
@pytest.mark.xfail
def test_remote_http():
    run(dpath("test_remote_http"))


@skip_on_windows
@connected
@pytest.mark.xfail
def test_remote_http_cluster():
    run(dpath("test_remote_http"), cluster=os.path.abspath(dpath("test14/qsub")))


def test_profile():
    run(dpath("test_profile"))


@skip_on_windows
@connected
def test_singularity():
    run(dpath("test_singularity"), use_singularity=True)


@skip_on_windows
def test_singularity_invalid():
    run(
        dpath("test_singularity"),
        targets=["invalid.txt"],
        use_singularity=True,
        shouldfail=True,
    )


@skip_on_windows
def test_singularity_module_invalid():
    run(
        dpath("test_singularity_module"),
        targets=["invalid.txt"],
        use_singularity=True,
        shouldfail=True,
    )


@skip_on_windows
@connected
def test_singularity_conda():
    run(
        dpath("test_singularity_conda"),
        use_singularity=True,
        use_conda=True,
        conda_frontend="conda",
    )


@skip_on_windows
@connected
def test_singularity_none():
    run(dpath("test_singularity_none"), use_singularity=True)


@skip_on_windows
@connected
def test_singularity_global():
    run(dpath("test_singularity_global"), use_singularity=True)


def test_issue612():
    run(dpath("test_issue612"), dryrun=True)


def test_bash():
    run(dpath("test_bash"))


def test_inoutput_is_path():
    run(dpath("test_inoutput_is_path"))


def test_archive():
    run(dpath("test_archive"), archive="workflow-archive.tar.gz")


def test_log_input():
    run(dpath("test_log_input"))


@skip_on_windows
@connected
def test_cwl():
    run(dpath("test_cwl"))


@skip_on_windows
@connected
def test_cwl_singularity():
    run(dpath("test_cwl"), use_singularity=True)


def test_issue805():
    run(dpath("test_issue805"), shouldfail=True)


def test_issue823_1():
    run(dpath("test_issue823_1"))


@skip_on_windows
def test_issue823_2():
    run(dpath("test_issue823_2"))


@skip_on_windows
def test_issue823_3():
    run(dpath("test_issue823_3"))


@skip_on_windows
def test_pathlib():
    run(dpath("test_pathlib"))


def test_pathlib_missing_file():
    run(dpath("test_pathlib_missing_file"), shouldfail=True)


@skip_on_windows
def test_group_jobs():
    run(dpath("test_group_jobs"), cluster="./qsub")


@skip_on_windows
def test_group_jobs_attempts():
    run(dpath("test_group_jobs_attempts"), cluster="./qsub", restart_times=2)


def assert_resources(resources: dict, **expected_resources):
    assert {res: resources[res] for res in expected_resources} == expected_resources


@skip_on_windows
def test_group_jobs_resources(mocker):
    spy = mocker.spy(GroupResources, "basic_layered")
    run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        cores=6,
        resources={"typo": 23, "mem_mb": 60000},
        group_components={0: 5},
        default_resources=DefaultResources(["mem_mb=0"]),
    )
    assert_resources(
        dict(spy.spy_return),
        _nodes=1,
        _cores=6,
        mem_mb=60000,
        runtime=420,
        fake_res=600,
        global_res=2000,
        disk_mb=2000,
    )


@skip_on_windows
def test_group_jobs_resources_with_max_threads(mocker):
    spy = mocker.spy(GroupResources, "basic_layered")
    run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        cores=6,
        resources={"mem_mb": 60000},
        max_threads=1,
        group_components={0: 5},
        default_resources=DefaultResources(["mem_mb=0"]),
    )
    assert_resources(
        dict(spy.spy_return),
        _nodes=1,
        _cores=5,
        mem_mb=60000,
        runtime=380,
        fake_res=1200,
        global_res=3000,
        disk_mb=3000,
    )


@skip_on_windows
def test_group_jobs_resources_with_limited_resources(mocker):
    spy = mocker.spy(GroupResources, "basic_layered")
    run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        cores=5,
        resources={"mem_mb": 10000},
        max_threads=1,
        group_components={0: 5},
        default_resources=DefaultResources(["mem_mb=0"]),
    )
    assert_resources(
        dict(spy.spy_return),
        _nodes=1,
        _cores=1,
        mem_mb=10000,
        runtime=700,
        fake_res=400,
        global_res=1000,
        disk_mb=1000,
    )


@skip_on_windows
def test_global_resource_limits_limit_scheduling_of_groups():
    # Note that in the snakefile, mem_mb is set as global (breaking the default) and
    # fake_res is set as local
    tmp = run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        cluster_status="./status_failed",
        cores=6,
        nodes=5,
        cleanup=False,
        resources={"typo": 23, "mem_mb": 50000, "fake_res": 200},
        group_components={0: 5, 1: 5},
        overwrite_groups={"a": 0, "a_1": 1, "b": 2, "c": 2},
        default_resources=DefaultResources(["mem_mb=0"]),
        shouldfail=True,
    )
    with (Path(tmp) / "qsub.log").open("r") as f:
        lines = [l for l in f.readlines() if not l == "\n"]
    assert len(lines) == 1
    shutil.rmtree(tmp)


@skip_on_windows
def test_new_resources_can_be_defined_as_local():
    # Test only works if both mem_mb and global_res are overwritten as local
    tmp = run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        cluster_status="./status_failed",
        cores=6,
        nodes=5,
        cleanup=False,
        resources={"typo": 23, "mem_mb": 50000, "fake_res": 200, "global_res": 1000},
        overwrite_resource_scopes={"mem_mb": "local", "global_res": "local"},
        group_components={0: 5, 1: 5},
        overwrite_groups={"a": 0, "a_1": 1, "b": 2, "c": 2},
        default_resources=DefaultResources(["mem_mb=0"]),
        shouldfail=True,
    )
    with (Path(tmp) / "qsub.log").open("r") as f:
        lines = [l for l in f.readlines() if not l == "\n"]
    assert len(lines) == 2
    shutil.rmtree(tmp)


@skip_on_windows
def test_resources_can_be_overwritten_as_global():
    # Test only works if fake_res overwritten as global
    tmp = run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        cluster_status="./status_failed",
        cores=6,
        nodes=5,
        cleanup=False,
        resources={"typo": 23, "fake_res": 200},
        overwrite_resource_scopes={"fake_res": "global"},
        group_components={0: 5, 1: 5},
        overwrite_groups={"a": 0, "a_1": 1, "b": 2, "c": 2},
        default_resources=DefaultResources(["mem_mb=0"]),
        shouldfail=True,
    )
    with (Path(tmp) / "qsub.log").open("r") as f:
        lines = [l for l in f.readlines() if not l == "\n"]
    assert len(lines) == 1
    shutil.rmtree(tmp)


@skip_on_windows
def test_scopes_submitted_to_cluster(mocker):
    from snakemake.executors import AbstractExecutor

    spy = mocker.spy(AbstractExecutor, "get_resource_scopes_args")
    run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        # cluster_status="./status_failed",
        overwrite_resource_scopes={"fake_res": "local"},
        max_threads=1,
        default_resources=DefaultResources(["mem_mb=0"]),
    )

    assert spy.spy_return == "--set-resource-scopes 'fake_res=local'"


@skip_on_windows
def test_resources_submitted_to_cluster(mocker):
    from snakemake.executors import AbstractExecutor

    spy = mocker.spy(AbstractExecutor, "get_resource_declarations_dict")
    run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        cores=6,
        resources={"mem_mb": 60000},
        max_threads=1,
        group_components={0: 5},
        default_resources=DefaultResources(["mem_mb=0"]),
    )

    assert_resources(
        spy.spy_return, mem_mb=60000, fake_res=1200, global_res=3000, disk_mb=3000
    )


@skip_on_windows
def test_excluded_resources_not_submitted_to_cluster(mocker):
    from snakemake.executors import AbstractExecutor

    spy = mocker.spy(AbstractExecutor, "get_resource_declarations_dict")
    run(
        dpath("test_group_jobs_resources"),
        cluster="./qsub",
        cores=6,
        resources={"mem_mb": 60000},
        max_threads=1,
        overwrite_resource_scopes={"fake_res": "excluded"},
        group_components={0: 5},
        default_resources=DefaultResources(["mem_mb=0"]),
    )
    assert_resources(spy.spy_return, mem_mb=60000, global_res=3000, disk_mb=3000)


@skip_on_windows
def test_group_job_resources_with_pipe(mocker):
    import copy
    from snakemake.executors import RealExecutor

    spy = mocker.spy(GroupResources, "basic_layered")

    # Cluster jobs normally get submitted with cores=all, but for testing purposes,
    # the system may not actually have enough cores to run the snakemake subprocess
    # (e.g. gh actions has only 2 cores). So cheat by patching over get_job_args
    get_job_args = copy.copy(RealExecutor.get_job_args)
    RealExecutor.get_job_args = lambda *args, **kwargs: get_job_args(
        *args, **{**kwargs, "cores": 6}
    )
    run(
        dpath("test_group_with_pipe"),
        cluster="./qsub",
        cores=6,
        resources={"mem_mb": 60000},
        group_components={0: 5},
        default_resources=DefaultResources(["mem_mb=0"]),
    )

    # revert to original version
    RealExecutor.get_job_args = get_job_args

    assert_resources(
        dict(spy.spy_return),
        _nodes=1,
        _cores=6,
        runtime=240,
        mem_mb=50000,
        disk_mb=1000,
    )


@skip_on_windows
def test_group_job_resources_with_pipe_with_too_much_constraint():
    run(
        dpath("test_group_with_pipe"),
        cluster="./qsub",
        cores=6,
        resources={"mem_mb": 20000},
        group_components={0: 5},
        shouldfail=True,
        default_resources=DefaultResources(["mem_mb=0"]),
    )


@skip_on_windows
def test_multicomp_group_jobs():
    run(
        dpath("test_multicomp_group_jobs"),
        cluster="./qsub",
        overwrite_groups={"a": "group0", "b": "group0"},
        group_components={"group0": 2},
    )


@skip_on_windows
def test_group_job_fail():
    run(dpath("test_group_job_fail"), cluster="./qsub", shouldfail=True)


@skip_on_windows  # Not supported, but could maybe be implemented. https://stackoverflow.com/questions/48542644/python-and-windows-named-pipes
def test_pipes():
    run(dpath("test_pipes"))


@skip_on_windows
def test_pipes_multiple():
    # see github issue #975
    run(dpath("test_pipes_multiple"), cores=5)


def test_pipes_fail():
    run(dpath("test_pipes_fail"), shouldfail=True)


def test_validate():
    run(dpath("test_validate"))


def test_validate_fail():
    run(
        dpath("test_validate"),
        configfiles=[dpath("test_validate/config.fail.yaml")],
        shouldfail=True,
    )


def test_issue854():
    # output and benchmark have inconsistent wildcards
    # this should fail when parsing
    run(dpath("test_issue854"), shouldfail=True)


@skip_on_windows
def test_issue850():
    run(dpath("test_issue850"), cluster="./qsub")


@skip_on_windows
def test_issue860():
    run(dpath("test_issue860"), cluster="./qsub", targets=["done"])


def test_issue894():
    run(dpath("test_issue894"))


def test_issue584():
    run(dpath("test_issue584"))


def test_issue912():
    run(dpath("test_issue912"))


@skip_on_windows
def test_job_properties():
    run(dpath("test_job_properties"), cluster="./qsub.py")


def test_issue916():
    run(dpath("test_issue916"))


@skip_on_windows
def test_issue930():
    run(dpath("test_issue930"), cluster="./qsub")


@skip_on_windows
def test_issue635():
    run(dpath("test_issue635"), use_conda=True, check_md5=False)


# TODO remove skip
@pytest.mark.skip(
    reason="Temporarily disable until the stable container image becomes available again."
)
@skip_on_windows
def test_convert_to_cwl():
    workdir = dpath("test_convert_to_cwl")
    # run(workdir, export_cwl=os.path.join(workdir, "workflow.cwl"))
    shell(
        "cd {workdir}; PYTHONPATH={src} python -m snakemake --export-cwl workflow.cwl",
        src=os.getcwd(),
    )
    shell("cd {workdir}; cwltool --singularity workflow.cwl")
    assert os.path.exists(os.path.join(workdir, "test.out"))


def test_issue1037():
    run(dpath("test_issue1037"), dryrun=True, cluster="qsub", targets=["Foo_A.done"])


def test_issue1046():
    run(dpath("test_issue1046"))


def test_checkpoints():
    run(dpath("test_checkpoints"))


def test_checkpoints_dir():
    run(dpath("test_checkpoints_dir"))


def test_issue1092():
    run(dpath("test_issue1092"))


@skip_on_windows
def test_issue1093():
    run(dpath("test_issue1093"), use_conda=True)


def test_issue958():
    run(dpath("test_issue958"), cluster="dummy", dryrun=True)


def test_issue471():
    run(dpath("test_issue471"))


def test_issue1085():
    run(dpath("test_issue1085"), shouldfail=True)


@skip_on_windows
def test_issue1083():
    run(dpath("test_issue1083"), use_singularity=True)


@skip_on_windows  # Fails with "The flag 'pipe' used in rule two is only valid for outputs
def test_pipes2():
    run(dpath("test_pipes2"))


def test_expand_flag():
    run(dpath("test_expand_flag"), shouldfail=True)


@skip_on_windows
def test_default_resources():
    from snakemake.resources import DefaultResources

    run(
        dpath("test_default_resources"),
        # use fractional defaults here to test whether they are correctly rounded
        default_resources=DefaultResources(
            ["mem_mb=max(2*input.size, 1000.1)", "disk_mb=max(2*input.size, 1000.2)"]
        ),
    )


@skip_on_windows  # TODO fix the windows case: it somehow does not consistently modify all temp env vars as desired
def test_tmpdir():
    # artificially set the tmpdir to an expected value
    run(dpath("test_tmpdir"), overwrite_resources={"a": {"tmpdir": "/tmp"}})


def test_tmpdir_default():
    # Do not check the content (OS and setup depdendent),
    # just check whether everything runs smoothly with the default.
    run(dpath("test_tmpdir"), check_md5=False)


def test_issue1284():
    run(dpath("test_issue1284"))


def test_issue1281():
    run(dpath("test_issue1281"))


def test_filegraph():
    workdir = dpath("test_filegraph")
    dot_path = os.path.join(workdir, "fg.dot")
    pdf_path = os.path.join(workdir, "fg.pdf")

    if ON_WINDOWS:
        shell.executable("bash")
        workdir = workdir.replace("\\", "/")
        dot_path = dot_path.replace("\\", "/")

    # make sure the calls work
    shell("cd {workdir};python -m snakemake --filegraph > {dot_path}")

    # make sure the output can be interpreted by dot
    with open(dot_path, "rb") as dot_file, open(pdf_path, "wb") as pdf_file:
        pdf_file.write(
            subprocess.check_output(["dot", "-Tpdf"], stdin=dot_file, cwd=workdir)
        )
    # make sure the generated pdf file is not empty
    assert os.stat(pdf_path).st_size > 0


def test_batch():
    from snakemake.dag import Batch

    run(dpath("test_batch"), batch=Batch("aggregate", 1, 2))


def test_batch_final():
    from snakemake.dag import Batch

    run(dpath("test_batch_final"), batch=Batch("aggregate", 1, 1))


def test_batch_fail():
    from snakemake.dag import Batch

    run(dpath("test_batch"), batch=Batch("aggregate", 2, 2), shouldfail=True)


def test_github_issue52():
    run(dpath("test_github_issue52"))
    run(dpath("test_github_issue52"), snakefile="other.smk")


@skip_on_windows
def test_github_issue78():
    run(dpath("test_github_issue78"), use_singularity=True)


def test_envvars():
    for var in ["TEST_ENV_VAR", "TEST_ENV_VAR2"]:
        try:
            del os.environ[var]
        except KeyError:
            pass
    run(dpath("test_envvars"), shouldfail=True)
    os.environ["TEST_ENV_VAR"] = "test"
    os.environ["TEST_ENV_VAR2"] = "test"
    run(dpath("test_envvars"), envvars=["TEST_ENV_VAR2"])


def test_github_issue105():
    run(dpath("test_github_issue105"))


def test_github_issue413():
    run(dpath("test_github_issue413"), no_tmpdir=True)


@skip_on_windows
def test_github_issue627():
    run(dpath("test_github_issue627"))


def test_github_issue727():
    run(dpath("test_github_issue727"))


@skip_on_windows
def test_github_issue988():
    run(dpath("test_github_issue988"))


@pytest.mark.skip(
    reason="ftp connections currently fail in github actions (TODO try again in the future)"
)
def test_github_issue1062():
    # old code failed in dry run
    run(dpath("test_github_issue1062"), dryrun=True)


def test_output_file_cache():
    test_path = dpath("test_output_file_cache")
    os.environ["SNAKEMAKE_OUTPUT_CACHE"] = "cache"
    run(test_path, cache=["a", "b"])
    run(test_path, cache=["invalid_multi"], targets="invalid1.txt", shouldfail=True)


@skip_on_windows
@pytest.mark.xfail(
    reason="moto currently fails with \"'_patch' object has no attribute 'is_local'\""
)
def test_output_file_cache_remote():
    test_path = dpath("test_output_file_cache_remote")
    os.environ["SNAKEMAKE_OUTPUT_CACHE"] = "cache"
    run(
        test_path,
        cache=["a", "b", "c"],
        default_remote_provider="S3Mocked",
        default_remote_prefix="test-remote-bucket",
    )


@connected
@zenodo
def test_remote_zenodo():
    run(dpath("test_remote_zenodo"))


def test_multiext():
    run(dpath("test_multiext"))


def test_core_dependent_threads():
    run(dpath("test_core_dependent_threads"))


@skip_on_windows
def test_env_modules():
    run(dpath("test_env_modules"), use_env_modules=True)


class TestParseCoresJobs:
    def run_test(self, func, ref):
        if ref is None:
            with pytest.raises(CliException):
                func()
            return
        assert func() == ref

    @pytest.mark.parametrize(
        ("input", "output"),
        [
            [(1, 1), (1, 1)],
            [(4, 4), (4, 4)],
            [(None, None), (1, 1)],
            [("all", "unlimited"), (available_cpu_count(), sys.maxsize)],
        ],
    )
    def test_no_exec(self, input, output):
        self.run_test(lambda: parse_cores_jobs(*input, True, False, False), output)
        # Test dryrun seperately
        self.run_test(lambda: parse_cores_jobs(*input, False, False, True), output)

    @pytest.mark.parametrize(
        ("input", "output"),
        [
            [(1, 1), (1, 1)],
            [(4, 4), (4, 4)],
            [(None, 1), (None, 1)],
            [(None, None), None],
            [(1, None), None],
            [("all", "unlimited"), (available_cpu_count(), sys.maxsize)],
        ],
    )
    def test_non_local_job(self, input, output):
        self.run_test(lambda: parse_cores_jobs(*input, False, True, False), output)

    @pytest.mark.parametrize(
        ("input", "output"),
        [
            [(1, 1), (1, None)],
            [(4, 4), (4, None)],
            [(None, 1), (1, None)],
            [(None, None), None],
            [(1, None), (1, None)],
            [(None, "all"), (available_cpu_count(), None)],
            [(None, "unlimited"), None],
            [("all", "unlimited"), (available_cpu_count(), None)],
        ],
    )
    def test_local_job(self, input, output):
        self.run_test(lambda: parse_cores_jobs(*input, False, False, False), output)


@skip_on_windows
@connected
def test_container():
    run(dpath("test_container"), use_singularity=True)


def test_linting():
    snakemake(
        snakefile=os.path.join(dpath("test14"), "Snakefile.nonstandard"), lint=True
    )


@skip_on_windows
def test_string_resources():
    from snakemake.resources import DefaultResources

    run(
        dpath("test_string_resources"),
        default_resources=DefaultResources(["gpu_model='nvidia-tesla-1000'"]),
        cluster="./qsub.py",
    )


def test_jupyter_notebook():
    run(dpath("test_jupyter_notebook"), use_conda=True)


def test_jupyter_notebook_draft():
    from snakemake.notebook import EditMode

    run(
        dpath("test_jupyter_notebook_draft"),
        use_conda=True,
        edit_notebook=EditMode(draft_only=True),
        targets=["results/result_intermediate.txt"],
        check_md5=False,
    )


def test_github_issue456():
    run(dpath("test_github_issue456"))


def test_scatter_gather():
    run(dpath("test_scatter_gather"), overwrite_scatter={"split": 2})


# SLURM tests go here, after successfull tests


def test_parsing_terminal_comment_following_statement():
    run(dpath("test_parsing_terminal_comment_following_statement"))


@skip_on_windows
def test_github_issue640():
    run(
        dpath("test_github_issue640"),
        targets=["Output/FileWithRights"],
        dryrun=True,
        cleanup=False,
    )


@skip_on_windows  # TODO check whether this might be enabled later
def test_generate_unit_tests():
    with tempfile.NamedTemporaryFile() as tmpfile:
        os.environ["UNIT_TEST_TMPFILE"] = tmpfile.name
        tmpdir = run(
            dpath("test_generate_unit_tests"),
            generate_unit_tests=".tests/unit",
            check_md5=False,
            cleanup=False,
        )
        sp.check_call(["pytest", ".tests", "-vs"], cwd=tmpdir)


@skip_on_windows
def test_metadata_migration():
    outpath = Path(
        "tests/test_metadata_migration/some/veryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryverylong"
    )
    os.makedirs(outpath, exist_ok=True)
    metapath = Path(
        "tests/test_metadata_migration/.snakemake/metadata/@c29tZS92ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5L3Zlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5L3Zlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcn/@l2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeS92ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnkvdmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnkvdmVyeXZlcnl2ZXJ5dmVy"
    )
    os.makedirs(metapath, exist_ok=True)
    exppath = Path(
        "tests/test_metadata_migration/expected-results/some/veryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryvery/veryveryveryveryveryveryveryveryveryveryveryveryverylong"
    )
    os.makedirs(exppath, exist_ok=True)
    with open(outpath / "path.txt", "w"):
        # generate empty file
        pass
    # generate artificial incomplete metadata in v1 format for migration
    with open(
        metapath
        / "eXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeWxvbmcvcGF0aC50eHQ=",
        "w",
    ) as meta:
        print('{"incomplete": true, "external_jobid": null}', file=meta)
    with open(exppath / "path.txt", "w") as out:
        print("updated", file=out)

    # run workflow, incomplete v1 metadata should be migrated and trigger rerun of the rule,
    # which will save different data than the output contained in the git repo.
    run(dpath("test_metadata_migration"), force_incomplete=True)


def test_paramspace():
    run(dpath("test_paramspace"))


def test_paramspace_single_wildcard():
    run(dpath("test_paramspace_single_wildcard"))


def test_github_issue806():
    run(dpath("test_github_issue806"), config=dict(src_lang="es", trg_lang="en"))


@skip_on_windows
def test_containerized():
    run(dpath("test_containerized"), use_conda=True, use_singularity=True)


@skip_on_windows
def test_containerize():
    run(dpath("test_conda"), containerize=True, check_results=False)


def test_long_shell():
    run(dpath("test_long_shell"))


def test_modules_all():
    run(dpath("test_modules_all"), targets=["a"])


def test_modules_all_exclude_1():
    # Fail due to conflicting rules
    run(dpath("test_modules_all_exclude"), shouldfail=True)


def test_modules_all_exclude_2():
    # Successed since the conflicting rules have been excluded
    run(
        dpath("test_modules_all_exclude"),
        snakefile="Snakefile_exclude",
        shouldfail=False,
    )


@skip_on_windows
def test_modules_prefix():
    run(dpath("test_modules_prefix"), targets=["a"])


@skip_on_windows
def test_modules_peppy():
    run(dpath("test_modules_peppy"), targets=["a"])


def test_modules_specific():
    run(dpath("test_modules_specific"), targets=["test_a"])


@skip_on_windows  # works in principle but the test framework modifies the target path separator
def test_modules_meta_wrapper():
    run(dpath("test_modules_meta_wrapper"), targets=["mapped/a.bam.bai"], dryrun=True)


def test_use_rule_same_module():
    run(dpath("test_use_rule_same_module"), targets=["test.out", "test2.out"])


def test_module_complex():
    run(dpath("test_module_complex"), dryrun=True)


def test_module_complex2():
    run(dpath("test_module_complex2"), dryrun=True)


@skip_on_windows
def test_module_no_prefixing_modified_paths():
    run(
        dpath("test_module_no_prefixing_modified_paths"),
        targets=["module2/test_final.txt"],
    )


def test_module_with_script():
    run(dpath("test_module_with_script"))


def test_module_worfklow_namespacing():
    run(dpath("test_module_workflow_snakefile_usage"))


@skip_on_windows  # No conda-forge version of pygraphviz for windows
def test_module_report():
    run(
        dpath("test_module_report"),
        report="report.html",
        report_stylesheet="custom-stylesheet.css",
        check_md5=False,
    )


def test_handover():
    run(dpath("test_handover"), resources={"mem_mb": 20})


@skip_on_windows  # test shell command not properly working
def test_source_path():
    run(dpath("test_source_path"), snakefile="workflow/Snakefile")


@only_on_windows
def test_filesep_windows_targets():
    run(
        dpath("test_filesep_windows"),
        targets=["subfolder/test2.out2", "subfolder/test1.out2"],
    )


@only_on_windows
def test_filesep_on_windows():
    run(dpath("test_filesep_windows"))


def test_set_resources():
    run(dpath("test_set_resources"), overwrite_resources={"a": {"a": 1, "b": "foo"}})


def test_github_issue1069():
    run(
        dpath("test_github_issue1069"),
        shellcmd="snakemake -c1 --resources mem_mb=16423",
    )


def test_touch_pipeline_with_temp_dir():
    # Issue #1028
    run(dpath("test_touch_pipeline_with_temp_dir"), forceall=True, touch=True)


def test_all_temp():
    run(dpath("test_all_temp"), all_temp=True)


def test_strict_mode():
    run(dpath("test_strict_mode"), shouldfail=True)


@needs_strace
def test_github_issue1158():
    run(dpath("test_github_issue1158"), cluster="./qsub.py")


def test_converting_path_for_r_script():
    run(dpath("test_converting_path_for_r_script"), cores=1)


def test_ancient_dag():
    run(dpath("test_ancient_dag"))


@skip_on_windows
def test_checkpoint_allowed_rules():
    run(dpath("test_checkpoint_allowed_rules"), targets=["c"], cluster="./qsub")


@skip_on_windows
def test_modules_ruledeps_inheritance():
    run(dpath("test_modules_ruledeps_inheritance"))


@skip_on_windows
def test_issue1331():
    # not guaranteed to fail, so let's try multiple times
    for i in range(10):
        run(dpath("test_issue1331"), cores=4)


@skip_on_windows
def test_conda_named():
    run(dpath("test_conda_named"), use_conda=True)


@skip_on_windows
def test_conda_function():
    run(dpath("test_conda_function"), use_conda=True, cores=1)


@skip_on_windows
def test_default_target():
    run(dpath("test_default_target"))


def test_cache_multioutput():
    run(dpath("test_cache_multioutput"), shouldfail=True)


@skip_on_windows
def test_github_issue1384():
    try:
        tmpdir = run(dpath("test_github_issue1384"), cleanup=False)
        shell(
            """
            cd {tmpdir}
            python -m snakemake --generate-unit-tests
            pytest -v .tests/unit
            """
        )
    finally:
        shutil.rmtree(tmpdir)


@skip_on_windows
def test_peppy():
    run(dpath("test_peppy"))


def test_template_engine():
    run(dpath("test_template_engine"))


def test_groupid_expand_local():
    run(dpath("test_groupid_expand"))


@skip_on_windows
def test_groupid_expand_cluster():
    run(dpath("test_groupid_expand_cluster"), cluster="./qsub", nodes=3)


@skip_on_windows
def test_service_jobs():
    run(dpath("test_service_jobs"), check_md5=False)


def test_incomplete_params():
    run(dpath("test_incomplete_params"), dryrun=True, printshellcmds=True)


@skip_on_windows
def test_github_issue261():
    run(dpath("test_github_issue261"), targets=["test1/target1/config1.done"])


@skip_on_windows  # no pipe support on windows
def test_pipe_depend():
    run(dpath("test_pipe_depend"), shouldfail=True)


@skip_on_windows  # no pipe support on windows
def test_pipe_depend_target_file():
    run(dpath("test_pipe_depend"), targets=["test.txt"], shouldfail=True)


@skip_on_windows  # platform independent issue
def test_github_issue1500():
    run(dpath("test_github_issue1500"), dryrun=True)


def test_github_issue1542():
    run(dpath("test_github_issue1542"), dryrun=True)


def test_github_issue1550():
    from snakemake.resources import DefaultResources

    run(
        dpath("test_github_issue1550"),
        resources={"mem_mb": 4000},
        default_resources=DefaultResources(
            ["mem_mb=max(2*input.size, 1000)", "disk_mb=max(2*input.size, 1000)"]
        ),
    )


def test_github_issue1498():
    run(dpath("test_github_issue1498"))


def test_cleanup_metadata_fail():
    run(dpath("test09"), cleanup_metadata=["xyz"])


@skip_on_windows  # same on win, no need to test
def test_github_issue1389():
    run(dpath("test_github_issue1389"), resources={"foo": 4}, shouldfail=True)


def test_github_issue2142():
    run(dpath("test_github_issue2142"))


def test_ensure_nonempty_fail():
    run(dpath("test_ensure"), targets=["a"], shouldfail=True)


def test_ensure_success():
    run(dpath("test_ensure"), targets=["b", "c"])


def test_ensure_checksum_fail():
    run(dpath("test_ensure"), targets=["d"], shouldfail=True)


@skip_on_windows
def test_github_issue1261():
    run(dpath("test_github_issue1261"), shouldfail=True, check_results=True)


def test_rule_inheritance_globals():
    run(
        dpath("test_rule_inheritance_globals"),
        report="report.html",
        targets=["foo.txt"],
        check_md5=False,
    )


def test_retries():
    run(dpath("test_retries"))


def test_retries_not_overriden():
    run(dpath("test_retries_not_overriden"), restart_times=3, shouldfail=True)


@skip_on_windows  # OS agnostic
def test_module_input_func():
    run(dpath("test_module_input_func"))


@skip_on_windows  # the testcase only has a linux-64 pin file
def test_conda_pin_file():
    run(dpath("test_conda_pin_file"), use_conda=True)


@skip_on_windows  # sufficient to test this on linux
def test_github_issue1618():
    run(dpath("test_github_issue1618"), cores=5)


def test_conda_python_script():
    run(dpath("test_conda_python_script"), use_conda=True)


@skip_on_windows
def test_github_issue1818():
    run(dpath("test_github_issue1818"), rerun_triggers="input")


@skip_on_windows  # not platform dependent
def test_match_by_wildcard_names():
    run(dpath("test_match_by_wildcard_names"))


@skip_on_windows  # not platform dependent
def test_github_issue929():
    # Huge thanks to Elmar Pruesse for providing this test case
    # and pointing to the problem in the code!
    # Huge thanks to Ronald Lehnigk for pointing me to the issue!
    run(dpath("test_github_issue929"), targets=["childrule_2"])


def test_github_issue1882():
    try:
        tmpdir = run(dpath("test_github_issue1882"), cleanup=False)
        run(tmpdir, forceall=True)
    finally:
        shutil.rmtree(tmpdir)


@skip_on_windows  # not platform dependent
def test_inferred_resources():
    run(dpath("test_inferred_resources"))


@skip_on_windows
def test_localrule():
    run(dpath("test_localrule"), targets=["1.txt", "2.txt"])


@skip_on_windows
def test_module_wildcard_constraints():
    run(dpath("test_module_wildcard_constraints"))
