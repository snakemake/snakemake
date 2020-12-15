__authors__ = ["Tobias Marschall", "Marcel Martin", "Johannes Köster"]
__copyright__ = "Copyright 2015-2020, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import sys
import uuid
import subprocess as sp
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

from common import *
from .conftest import skip_on_windows, ON_WINDOWS


def test_list_untracked():
    run(dpath("test_list_untracked"))


xfail_permissionerror_on_win = (
    pytest.mark.xfail(raises=PermissionError) if ON_WINDOWS else lambda x: x
)

# Must fail on Windows with PermissionError since the tempfile.TemporaryDirectory
# can't clean up the protected files generated in the test
@xfail_permissionerror_on_win
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
    run(dpath("test14"), snakefile="Snakefile.nonstandard", cluster="./qsub")


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
        targets=["downstream", "symlinked_input", "child_to_input"],
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


@skip_on_windows  # No conda-forge version of pygraphviz for windows
def test_report():
    run(
        dpath("test_report"),
        report="report.html",
        report_stylesheet="custom-stylesheet.css",
        check_md5=False,
    )


@skip_on_windows  # No conda-forge version of pygraphviz for windows
def test_report_zip():
    run(dpath("test_report_zip"), report="report.zip", check_md5=False)


@skip_on_windows  # No conda-forge version of pygraphviz for windows
def test_report_dir():
    run(dpath("test_report_dir"), report="report.zip", check_md5=False)


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


def test_config():
    run(dpath("test_config"))


def test_update_config():
    run(dpath("test_update_config"))


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


@skip_on_windows  # Symbolic link privileges needed to work
def test_shadow_prefix():
    run(dpath("test_shadow_prefix"), shadow_prefix="shadowdir")


@skip_on_windows
def test_shadow_prefix_qsub():
    run(dpath("test_shadow_prefix"), shadow_prefix="shadowdir", cluster="./qsub")


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


@skip_on_windows  # test uses bwa which is non windows
def test_conda():
    run(dpath("test_conda"), use_conda=True)


@skip_on_windows  # test uses bwa which is non windows
def test_upstream_conda():
    run(dpath("test_conda"), use_conda=True, conda_frontend="conda")


@skip_on_windows  # Conda support is partly broken on Win
def test_conda_custom_prefix():
    run(
        dpath("test_conda_custom_prefix"),
        use_conda=True,
        conda_prefix="custom",
        set_pythonpath=False,
    )


@skip_on_windows  # Conda support is partly broken on Win
def test_wrapper():
    run(dpath("test_wrapper"), use_conda=True)


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

    The qsub in the sub directory will fail the first time and succeed the
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
        User project to bill for download. None will not provide project for
        requester-pays as is the usual default
    requesting_url: str
        URL of bucket to download. Default will match expected output, but is a
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
@connected
def test_singularity_conda():
    run(
        dpath("test_singularity_conda"),
        use_singularity=True,
        use_conda=True,
        conda_frontend="conda",
    )


def test_issue612():
    run(dpath("test_issue612"), dryrun=True)


def test_bash():
    run(dpath("test_bash"))


def test_inoutput_is_path():
    run(dpath("test_inoutput_is_path"))


# Fails on Windows with PermissionError when test system tries to
# clean the conda environment
@xfail_permissionerror_on_win
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


@skip_on_windows
def test_pathlib():
    run(dpath("test_pathlib"))


def test_pathlib_missing_file():
    run(dpath("test_pathlib_missing_file"), shouldfail=True)


@skip_on_windows
def test_group_jobs():
    run(dpath("test_group_jobs"), cluster="./qsub")


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
        default_resources=DefaultResources(
            ["mem_mb=max(2*input.size, 1000)", "disk_mb=max(2*input.size, 1000)"]
        ),
    )


def test_issue1284():
    run(dpath("test_issue1284"))


def test_issue1281():
    run(dpath("test_issue1281"))


@skip_on_windows  # Currently no workable pygraphviz package
def test_filegraph():
    workdir = dpath("test_filegraph")
    dot_path = os.path.abspath("fg.dot")
    pdf_path = "fg.pdf"

    # make sure the calls work
    shell("cd {workdir}; python -m snakemake --filegraph > {dot_path}")

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
    run(dpath("test_envvars"), shouldfail=True)
    os.environ["TEST_ENV_VAR"] = "test"
    os.environ["TEST_ENV_VAR2"] = "test"
    run(dpath("test_envvars"), envvars=["TEST_ENV_VAR2"])


def test_github_issue105():
    run(dpath("test_github_issue105"))


def test_output_file_cache():
    test_path = dpath("test_output_file_cache")
    os.environ["SNAKEMAKE_OUTPUT_CACHE"] = os.path.join(test_path, "cache")
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


def test_multiext():
    run(dpath("test_multiext"))


def test_core_dependent_threads():
    run(dpath("test_core_dependent_threads"))


@skip_on_windows
def test_env_modules():
    run(dpath("test_env_modules"), use_env_modules=True)


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


@skip_on_windows  # currently fails on windows. Please help fix.
def test_jupyter_notebook():
    run(dpath("test_jupyter_notebook"), use_conda=True)


def test_github_issue456():
    run(dpath("test_github_issue456"))


def test_scatter_gather():
    run(dpath("test_scatter_gather"), overwrite_scatter={"split": 2})


@skip_on_windows
def test_github_issue640():
    run(
        dpath("test_github_issue640"),
        targets=["Output/FileWithRights"],
        dryrun=True,
        cleanup=False,
    )


def test_generate_unit_tests():
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
        / "eXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeXZlcnl2ZXJ5dmVyeWxvbmcvcGF0aC50eHQ\=",
        "w",
    ) as meta:
        print('{"incomplete": true, "external_jobid": null}', file=meta)
    with open(exppath / "path.txt", "w") as out:
        print("updated", file=out)

    # run workflow, incomplete v1 metadata should be migrated and trigger rerun of the rule,
    # which will save different data than the output contained in the git repo.
    run(dpath("test_metadata_migration"), force_incomplete=True)
