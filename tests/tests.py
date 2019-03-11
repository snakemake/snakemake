__authors__ = ["Tobias Marschall", "Marcel Martin", "Johannes Köster"]
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import sys
import os
import shutil
from os.path import join
from subprocess import call
import tempfile
import hashlib
import urllib
from shutil import rmtree, which
from shlex import quote
import pytest
import subprocess

from snakemake import snakemake
from snakemake.shell import shell


if not which("snakemake"):
    raise Exception("snakemake not in PATH. For testing, install snakemake with "
                    "'pip install -e .'. You should do this in a separate environment "
                    "(via conda or virtualenv).")


def dpath(path):
    """get path to a data file (relative to the directory this
	test lives in)"""
    return os.path.realpath(join(os.path.dirname(__file__), path))


def md5sum(filename):
    data = open(filename, 'rb').read()
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
    return "GCLOUD_SERVICE_KEY" in os.environ


def has_gcloud_cluster():
    return "GCLOUD_CLUSTER" in os.environ

gcloud = pytest.mark.skipif(not is_connected()
                                 or not has_gcloud_service_key()
                                 or not has_gcloud_cluster(),
                                 reason="Skipping GCLOUD tests because not on "
                                        "CI, no inet connection or not logged "
                                        "in to gcloud.")

connected = pytest.mark.skipif(not is_connected(), reason="no internet connection")

ci = pytest.mark.skipif(not is_ci(), reason="not in CI")

def copy(src, dst):
    if os.path.isdir(src):
        shutil.copytree(src, os.path.join(dst, os.path.basename(src)))
    else:
        shutil.copy(src, dst)


def run(path,
        shouldfail=False,
        snakefile="Snakefile",
        subpath=None,
        no_tmpdir=False,
        check_md5=True, cores=3, **params):
    """
    Test the Snakefile in path.
    There must be a Snakefile in the path and a subdirectory named
    expected-results.
    """
    results_dir = join(path, 'expected-results')
    snakefile = join(path, snakefile)
    assert os.path.exists(snakefile)
    assert os.path.exists(results_dir) and os.path.isdir(
        results_dir), '{} does not exist'.format(results_dir)
    with tempfile.TemporaryDirectory(prefix="snakemake-") as tmpdir:
        config = {}
        # handle subworkflow
        if subpath is not None:
            # set up a working directory for the subworkflow and pass it in `config`
            # for now, only one subworkflow is supported
            assert os.path.exists(subpath) and os.path.isdir(
                subpath), '{} does not exist'.format(subpath)
            subworkdir = os.path.join(tmpdir, "subworkdir")
            os.mkdir(subworkdir)
            # copy files
            for f in os.listdir(subpath):
                copy(os.path.join(subpath, f), subworkdir)
            config['subworkdir'] = subworkdir

        # copy files
        for f in os.listdir(path):
            print(f)
            copy(os.path.join(path, f), tmpdir)

        # run snakemake
        success = snakemake(snakefile,
                            cores=cores,
                            workdir=path if no_tmpdir else tmpdir,
                            stats="stats.txt",
                            config=config, **params)
        if shouldfail:
            assert not success, "expected error on execution"
        else:
            assert success, "expected successful execution"
            for resultfile in os.listdir(results_dir):
                if resultfile == ".gitignore" or not os.path.isfile(
                    os.path.join(results_dir, resultfile)):
                    # this means tests cannot use directories as output files
                    continue
                targetfile = join(tmpdir, resultfile)
                expectedfile = join(results_dir, resultfile)
                assert os.path.exists(
                    targetfile), 'expected file "{}" not produced'.format(
                        resultfile)
                if check_md5:
                    # if md5sum(targetfile) != md5sum(expectedfile):
                    #     import pdb; pdb.set_trace()
                    assert md5sum(targetfile) == md5sum(
                        expectedfile), 'wrong result produced for file "{}"'.format(
                            resultfile)

def test_delete_all_output():
    run(dpath("test_delete_all_output"))

def test_issue956():
    run(dpath("test_issue956"))

def test01():
    run(dpath("test01"))


def test02():
    run(dpath("test02"))


def test03():
    run(dpath("test03"), targets=['test.out'])


def test04():
    run(dpath("test04"), targets=['test.out'])


def test05():
    run(dpath("test05"))


def test06():
    run(dpath("test06"), targets=['test.bla.out'])


def test07():
    run(dpath("test07"), targets=['test.out', 'test2.out'])


def test08():
    run(dpath("test08"), targets=['test.out', 'test2.out'])


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


def test14():
    run(dpath("test14"), snakefile="Snakefile.nonstandard", cluster="./qsub")


def test15():
    run(dpath("test15"))

def test_directory():
    run(dpath("test_directory"), targets=['downstream', 'symlinked_input', "child_to_input", "some/dir-child", "some/shadow"])
    run(dpath("test_directory"), targets=['file_expecting_dir'], shouldfail = True)
    run(dpath("test_directory"), targets=['dir_expecting_file'], shouldfail = True)
    run(dpath("test_directory"), targets=['child_to_other'], shouldfail = True)

def test_ancient():
    run(dpath("test_ancient"), targets=['D', 'old_file'])

def test_list_untracked():
    run(dpath("test_list_untracked"))

def test_report():
    run(dpath("test_report"), report="report.html", check_md5=False)


def test_dynamic():
    run(dpath("test_dynamic"))


def test_params():
    run(dpath("test_params"))


def test_same_wildcard():
    run(dpath("test_same_wildcard"))


def test_conditional():
    run(dpath("test_conditional"),
        targets="test.out test.0.out test.1.out test.2.out".split())


def test_unpack_dict():
    run(dpath("test_unpack_dict"))


def test_unpack_list():
    run(dpath("test_unpack_list"))


def test_shell():
    run(dpath("test_shell"))


def test_temp():
    run(dpath("test_temp"),
        cluster="./qsub",
        targets="test.realigned.bam".split())


def test_keyword_list():
    run(dpath("test_keyword_list"))


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


def test_benchmark():
    run(dpath("test_benchmark"), check_md5=False)


def test_temp_expand():
    run(dpath("test_temp_expand"))


def test_wildcard_count_ambiguity():
    run(dpath("test_wildcard_count_ambiguity"))


def test_srcdir():
    run(dpath("test_srcdir"))


def test_multiple_includes():
    run(dpath("test_multiple_includes"))


def test_yaml_config():
    run(dpath("test_yaml_config"))


# TODO reenable once S3Mocked works with boto3
# def test_remote():
#     try:
#         import moto
#         import boto3
#         import filechunkio
#
#         # only run the remote file test if the dependencies
#         # are installed, otherwise do nothing
#         run(dpath("test_remote"), cores=1)
#     except ImportError:
#         pass


def test_cluster_sync():
    run(dpath("test14"),
        snakefile="Snakefile.nonstandard",
        cluster_sync="./qsub")

@pytest.mark.skip(reason="This does not work reliably in CircleCI.")
def test_symlink_temp():
    run(dpath("test_symlink_temp"), shouldfail=True)


def test_empty_include():
    run(dpath("test_empty_include"))


def test_script():
    run(dpath("test_script"))


def test_shadow():
    run(dpath("test_shadow"))

def test_shadow_prefix():
    run(dpath("test_shadow_prefix"), shadow_prefix = "shadowdir")
    run(dpath("test_shadow_prefix"), shadow_prefix = "shadowdir", cluster="./qsub")

def test_until():
    run(dpath("test_until"),
        until=["leveltwo_first", # rule name
               "leveltwo_second.txt", # file name
               "second_wildcard"]) # wildcard rule


def test_omitfrom():
    run(dpath("test_omitfrom"),
        omit_from=["leveltwo_first", # rule name
                   "leveltwo_second.txt", # file name
                   "second_wildcard"]) # wildcard rule


def test_nonstr_params():
    run(dpath("test_nonstr_params"))


def test_delete_output():
    run(dpath("test_delete_output"), cores=1)


def test_input_generator():
    run(dpath("test_input_generator"))


def test_symlink_time_handling():
    #See Snakefile for notes on why this fails on some systems
    if os.utime in os.supports_follow_symlinks:
        run(dpath("test_symlink_time_handling"))


def test_issue328():
    try:
        import pytools
        run(dpath("test_issue328"), forcerun=["split"])
    except ImportError:
        # skip test if import fails
        pass


def test_conda():
    if conda_available():
        run(dpath("test_conda"), use_conda=True)


def test_conda_custom_prefix():
    if conda_available():
        run(dpath("test_conda_custom_prefix"),
            use_conda=True, conda_prefix="custom")


def test_wrapper():
    if conda_available():
        run(dpath("test_wrapper"), use_conda=True)


def conda_available():
    return which("conda")


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
    run(dpath("test_spaces_in_fnames"),
        # cluster="./qsub",
        targets=["test bam file realigned.bam"],
        verbose=True,
        printshellcmds=True)


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
def test_remote_irods():
    if os.environ.get("CI") == "true":
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
    run(dpath("test (with parentheses)"))

def test_dup_out_patterns():
    """Duplicate output patterns should emit an error

    Duplicate output patterns can be detected on the rule level
    """
    run(dpath("test_dup_out_patterns"), shouldfail=True)


def test_restartable_job_cmd_exit_1():
    """Test the restartable job feature on ``exit 1``

    The shell snippet in the Snakemake file will fail the first time
    and succeed the second time.
    """
    # Even two consecutive times should fail as files are cleared
    run(dpath("test_restartable_job_cmd_exit_1"), cluster="./qsub",
        restart_times=0, shouldfail=True)
    run(dpath("test_restartable_job_cmd_exit_1"), cluster="./qsub",
        restart_times=0, shouldfail=True)
    # Restarting once is enough
    run(dpath("test_restartable_job_cmd_exit_1"), cluster="./qsub",
        restart_times=1, printshellcmds=True)


def test_restartable_job_qsub_exit_1():
    """Test the restartable job feature when qsub fails

    The qsub in the sub directory will fail the first time and succeed the
    second time.
    """
    # Even two consecutive times should fail as files are cleared
    run(dpath("test_restartable_job_qsub_exit_1"), cluster="./qsub",
        restart_times=0, shouldfail=True)
    run(dpath("test_restartable_job_qsub_exit_1"), cluster="./qsub",
        restart_times=0, shouldfail=True)
    # Restarting once is enough
    run(dpath("test_restartable_job_qsub_exit_1"), cluster="./qsub",
        restart_times=1, shouldfail=False)


def test_threads():
    run(dpath("test_threads"), cores=20)


def test_dynamic_temp():
    run(dpath("test_dynamic_temp"))


# TODO this currently hangs. Has to be investigated (issue #660).
#def test_ftp_immediate_close():
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

@pytest.mark.skip(reason="moto seems to be broken currently")
def test_default_remote():
     run(dpath("test_default_remote"),
         cores=1,
         default_remote_provider="S3Mocked",
         default_remote_prefix="test-remote-bucket")


def test_run_namedlist():
    run(dpath("test_run_namedlist"))


@connected
@ci
def test_remote_gs():
    if not "CI" in os.environ:
        run(dpath("test_remote_gs"))
    else:
        print("skipping test_remove_gs in CI")


@connected
def test_remote_log():
    run(dpath("test_remote_log"), shouldfail=True)


@connected
@pytest.mark.xfail
def test_remote_http():
    run(dpath("test_remote_http"))


@connected
@pytest.mark.xfail
def test_remote_http_cluster():
    run(dpath("test_remote_http"), cluster=os.path.abspath(dpath("test14/qsub")))

def test_profile():
    run(dpath("test_profile"))


@connected
def test_singularity():
    run(dpath("test_singularity"), use_singularity=True)

def test_singularity_invalid():
    run(dpath("test_singularity"), targets=["invalid.txt"], use_singularity=True, shouldfail=True)

@connected
def test_singularity_conda():
    run(dpath("test_singularity_conda"), use_singularity=True, use_conda=True)


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


@pytest.fixture(scope="module")
def gcloud_cluster():
    class Cluster:
        def __init__(self):
            self.cluster = os.environ["GCLOUD_CLUSTER"]
            self.bucket_name = 'snakemake-testing-{}'.format(self.cluster)

            shell("""
            $GCLOUD container clusters create {self.cluster} --num-nodes 3 --scopes storage-rw --zone us-central1-a --machine-type f1-micro
            $GCLOUD container clusters get-credentials {self.cluster} --zone us-central1-a
            $GSUTIL mb gs://{self.bucket_name}
            """)

        def delete(self):
            shell("""
            $GCLOUD container clusters delete {self.cluster} --zone us-central1-a --quiet || true
            $GSUTIL rm -r gs://{self.bucket_name} || true
            """)

        def run(self, test="test_kubernetes", **kwargs):
            try:
                run(dpath(test),
                    kubernetes="default",
                    default_remote_provider="GS",
                    default_remote_prefix=self.bucket_name,
                    no_tmpdir=True,
                    **kwargs)
            except Exception as e:
                shell("for p in `kubectl get pods | grep ^snakejob- | cut -f 1 -d ' '`; do kubectl logs $p; done")
                raise e

        def reset(self):
            shell('$GSUTIL rm -r gs://{self.bucket_name}/* || true')
    cluster = Cluster()
    yield cluster
    cluster.delete()



@gcloud
@pytest.mark.xfail
def test_gcloud_plain(gcloud_cluster):
    gcloud_cluster.reset()
    gcloud_cluster.run()


@gcloud
@pytest.mark.skip(reason="need a faster cloud compute instance to run this")
def test_gcloud_conda(gcloud_cluster):
    gcloud_cluster.reset()
    gcloud_cluster.run(use_conda=True)


@gcloud
@pytest.mark.skip(reason="need a faster cloud compute instance to run this")
def test_gcloud_singularity(gcloud_cluster):
    gcloud_cluster.reset()
    gcloud_cluster.run(use_singularity=True)


@gcloud
@pytest.mark.skip(reason="need a faster cloud compute instance to run this")
def test_gcloud_conda_singularity(gcloud_cluster):
    gcloud_cluster.reset()
    gcloud_cluster.run(use_singularity=True, use_conda=True)


@gcloud()
@pytest.mark.skip(reason="need a faster cloud compute instance to run this")
def test_issue1041(gcloud_cluster):
    gcloud_cluster.reset()
    gcloud_cluster.run(test="test_issue1041")


@connected
def test_cwl():
    run(dpath("test_cwl"))


@connected
def test_cwl_singularity():
    run(dpath("test_cwl"), use_singularity=True)


def test_issue805():
    run(dpath("test_issue805"), shouldfail=True)


def test_group_jobs():
    run(dpath("test_group_jobs"), cluster="./qsub")


def test_group_job_fail():
    run(dpath("test_group_job_fail"), cluster="./qsub", shouldfail=True)


def test_pipes():
    run(dpath("test_pipes"))


def test_pipes_fail():
    run(dpath("test_pipes_fail"), shouldfail=True)


def test_validate():
    run(dpath("test_validate"))


def test_validate_fail():
    run(dpath("test_validate"), configfile=dpath("test_validate/config.fail.yaml"), shouldfail=True)


def test_issue854():
    # output and benchmark have inconsistent wildcards
    # this should fail when parsing
    run(dpath("test_issue854"), shouldfail=True)


def test_issue850():
    run(dpath("test_issue850"), cluster="./qsub")


def test_issue860():
    run(dpath("test_issue860"), cluster="./qsub", targets=["done"])


def test_issue894():
    run(dpath("test_issue894"))

def test_issue584():
    run(dpath("test_issue584"))

def test_issue912():
    run(dpath("test_issue912"))

def test_job_properties():
    run(dpath("test_job_properties"), cluster="./qsub.py")

def test_issue916():
    run(dpath("test_issue916"))

def test_issue930():
    run(dpath("test_issue930"), cluster="./qsub")

def test_issue635():
    run(dpath("test_issue635"), use_conda=True, check_md5=False)

def test_convert_to_cwl():
    workdir = dpath("test_convert_to_cwl")
    #run(workdir, export_cwl=os.path.join(workdir, "workflow.cwl"))
    subprocess.check_call(["snakemake", "--export-cwl" , "workflow.cwl"], cwd=workdir)
    subprocess.check_call(["cwltool", "workflow.cwl"], cwd=workdir)
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

def test_issue1093():
    run(dpath("test_issue1093"), use_conda=True)

def test_issue958():
    run(dpath("test_issue958"), cluster="dummy", dryrun=True)

def test_issue471():
    run(dpath("test_issue471"))

def test_issue1085():
    run(dpath("test_issue1085"), shouldfail=True)

def test_issue1083():
    run(dpath("test_issue1083"), use_singularity=True)
