__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import abstractmethod
import os
import sys
import contextlib
import time
import datetime
import json
import textwrap
import stat
import shutil
import shlex
import threading
import concurrent.futures
import subprocess
import signal
import tempfile
import threading
from functools import partial
from itertools import chain
from collections import namedtuple
from snakemake.executors.common import format_cli_arg, format_cli_pos_arg, join_cli_args
from snakemake.io import _IOFile
import random
import base64
import uuid
import re
import math

from snakemake.jobs import Job
from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.stats import Stats
from snakemake.utils import format, Unformattable, makedirs
from snakemake.io import get_wildcard_names, Wildcards
from snakemake.exceptions import print_exception, get_exception_origin
from snakemake.exceptions import format_error, RuleException, log_verbose_traceback
from snakemake.exceptions import (
    ProtectedOutputException,
    WorkflowError,
    ImproperShadowException,
    SpawnedJobError,
    CacheMissException,
)
from snakemake.common import (
    Mode,
    __version__,
    get_container_image,
    get_uuid,
    lazy_property,
)


# TODO move each executor into a separate submodule


def sleep():
    # do not sleep on CI. In that case we just want to quickly test everything.
    if os.environ.get("CI") != "true":
        time.sleep(10)
    else:
        time.sleep(1)


class AbstractExecutor:
    def __init__(
        self,
        workflow,
        dag,
        printreason=False,
        quiet=False,
        printshellcmds=False,
        printthreads=True,
        keepincomplete=False,
    ):
        self.workflow = workflow
        self.dag = dag
        self.quiet = quiet
        self.printreason = printreason
        self.printshellcmds = printshellcmds
        self.printthreads = printthreads
        self.latency_wait = workflow.latency_wait
        self.keepincomplete = keepincomplete

    def get_default_remote_provider_args(self):
        return join_cli_args(
            [
                self.workflow_property_to_arg("default_remote_prefix"),
                self.workflow_property_to_arg("default_remote_provider", attr="name"),
            ]
        )

    def get_set_resources_args(self):
        return format_cli_arg(
            "--set-resources",
            [
                f"{rule}:{name}={value}"
                for rule, res in self.workflow.overwrite_resources.items()
                for name, value in res.items()
            ],
            skip=not self.workflow.overwrite_resources,
        )

    def get_default_resources_args(self, default_resources=None):
        default_resources = default_resources or self.workflow.default_resources
        return format_cli_arg("--default-resources", default_resources.args)

    def run_jobs(self, jobs, callback=None, submit_callback=None, error_callback=None):
        """Run a list of jobs that is ready at a given point in time.

        By default, this method just runs each job individually.
        This method can be overwritten to submit many jobs in a more efficient way than one-by-one.
        Note that in any case, for each job, the callback functions have to be called individually!
        """
        for job in jobs:
            self.run(
                job,
                callback=callback,
                submit_callback=submit_callback,
                error_callback=error_callback,
            )

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        """Run a specific job or group job."""
        self._run(job)
        callback(job)

    def shutdown(self):
        pass

    def cancel(self):
        pass

    def _run(self, job):
        job.check_protected_output()
        self.printjob(job)

    def rule_prefix(self, job):
        return "local " if job.is_local else ""

    def printjob(self, job):
        job.log_info(skip_dynamic=True)

    def print_job_error(self, job, msg=None, **kwargs):
        job.log_error(msg, **kwargs)

    def handle_job_success(self, job):
        pass

    def handle_job_error(self, job):
        pass


class DryrunExecutor(AbstractExecutor):
    def printjob(self, job):
        super().printjob(job)
        if job.is_group():
            for j in job.jobs:
                self.printcache(j)
        else:
            self.printcache(job)

    def printcache(self, job):
        if self.workflow.is_cached_rule(job.rule):
            if self.workflow.output_file_cache.exists(job):
                logger.info(
                    "Output file {} will be obtained from global between-workflow cache.".format(
                        job.output[0]
                    )
                )
            else:
                logger.info(
                    "Output file {} will be written to global between-workflow cache.".format(
                        job.output[0]
                    )
                )


class RealExecutor(AbstractExecutor):
    def __init__(
        self,
        workflow,
        dag,
        printreason=False,
        quiet=False,
        printshellcmds=False,
        assume_shared_fs=True,
        keepincomplete=False,
    ):
        super().__init__(
            workflow,
            dag,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            keepincomplete=keepincomplete,
        )
        self.assume_shared_fs = assume_shared_fs
        self.stats = Stats()
        self.snakefile = workflow.main_snakefile

    def register_job(self, job):
        job.register()

    def _run(self, job, callback=None, error_callback=None):
        super()._run(job)
        self.stats.report_job_start(job)

        try:
            self.register_job(job)
        except IOError as e:
            logger.info(
                "Failed to set marker file for job started ({}). "
                "Snakemake will work, but cannot ensure that output files "
                "are complete in case of a kill signal or power loss. "
                "Please ensure write permissions for the "
                "directory {}".format(e, self.workflow.persistence.path)
            )

    def handle_job_success(
        self,
        job,
        upload_remote=True,
        handle_log=True,
        handle_touch=True,
        ignore_missing_output=False,
    ):
        job.postprocess(
            upload_remote=upload_remote,
            handle_log=handle_log,
            handle_touch=handle_touch,
            ignore_missing_output=ignore_missing_output,
            latency_wait=self.latency_wait,
            assume_shared_fs=self.assume_shared_fs,
            keep_metadata=self.workflow.keep_metadata,
        )
        self.stats.report_job_end(job)

    def handle_job_error(self, job, upload_remote=True):
        job.postprocess(
            error=True,
            assume_shared_fs=self.assume_shared_fs,
            latency_wait=self.latency_wait,
        )

    def workflow_property_to_arg(
        self, property, flag=None, quote=True, skip=False, invert=False, attr=None
    ):
        if skip:
            return ""

        value = getattr(self.workflow, property)

        if value is not None and attr is not None:
            value = getattr(value, attr)

        if flag is None:
            flag = f"--{property.replace('_', '-')}"

        if invert and isinstance(value, bool):
            value = not value

        return format_cli_arg(flag, value, quote=quote)

    @lazy_property
    def general_args(self):
        """Return a string to add to self.exec_job that includes additional
        arguments from the command line. This is currently used in the
        ClusterExecutor and CPUExecutor, as both were using the same
        code. Both have base class of the RealExecutor.
        """
        w2a = self.workflow_property_to_arg

        return join_cli_args(
            [
                "--force",
                "--keep-target-files",
                "--keep-remote",
                "--max-inventory-time 0",
                "--nocolor",
                "--notemp",
                "--no-hooks",
                "--nolock",
                "--ignore-incomplete",
                w2a("cleanup_scripts", flag="--skip-script-cleanup"),
                w2a("shadow_prefix"),
                w2a("use_conda"),
                w2a("conda_frontend"),
                w2a("conda_prefix"),
                w2a("conda_base_path", skip=not self.assume_shared_fs),
                w2a("use_singularity"),
                w2a("singularity_prefix"),
                w2a("singularity_args"),
                w2a("execute_subworkflows", flag="--no-subworkflows", invert=True),
                w2a("max_threads"),
                w2a("use_env_modules", flag="--use-envmodules"),
                w2a("keep_metadata", flag="--drop-metadata", invert=True),
                w2a("wrapper_prefix"),
                w2a("overwrite_threads", flag="--set-threads"),
                w2a("overwrite_scatter", flag="--set-scatter"),
                w2a("local_groupid", skip=self.job_specific_local_groupid),
                w2a("conda_not_block_search_path_envvars"),
                w2a("overwrite_configfiles", flag="--configfiles"),
                w2a("config_args", flag="--config"),
                w2a("printshellcmds"),
                w2a("latency_wait"),
                w2a("scheduler_type", flag="--scheduler"),
                format_cli_arg(
                    "--scheduler-solver-path",
                    os.path.dirname(sys.executable),
                    skip=not self.assume_shared_fs,
                ),
                self.get_set_resources_args(),
                self.get_default_remote_provider_args(),
                self.get_default_resources_args(),
                self.get_workdir_arg(),
                format_cli_arg("--mode", self.get_exec_mode()),
            ]
        )

    def get_workdir_arg(self):
        return self.workflow_property_to_arg("overwrite_workdir", flag="--directory")

    def get_job_args(self, job, **kwargs):
        return join_cli_args(
            [
                format_cli_pos_arg(kwargs.get("target", self.get_job_targets(job))),
                # Restrict considered rules for faster DAG computation.
                # This does not work for updated jobs because they need
                # to be updated in the spawned process as well.
                format_cli_arg(
                    "--allowed-rules",
                    job.rules,
                    quote=False,
                    skip=job.is_branched or job.is_updated,
                ),
                # Ensure that a group uses its proper local groupid.
                format_cli_arg("--local-groupid", job.jobid, skip=not job.is_group()),
                format_cli_arg("--cores", kwargs.get("cores", self.cores)),
                format_cli_arg("--attempt", job.attempt),
                format_cli_arg("--force-use-threads", not job.is_group()),
            ]
        )

    @property
    def job_specific_local_groupid(self):
        return True

    def get_snakefile(self):
        return self.snakefile

    def get_job_targets(self, job):
        return job.get_targets()

    @abstractmethod
    def get_python_executable(self):
        ...

    @abstractmethod
    def get_exec_mode(self):
        ...

    @abstractmethod
    def get_envvar_declarations(self):
        ...

    def get_job_exec_prefix(self, job):
        return ""

    def get_job_exec_suffix(self, job):
        return ""

    def format_job_exec(self, job):
        prefix = self.get_job_exec_prefix(job)
        if prefix:
            prefix += " &&"
        suffix = self.get_job_exec_suffix(job)
        if suffix:
            suffix = f"&& {suffix}"
        return join_cli_args(
            [
                prefix,
                self.get_envvar_declarations(),
                self.get_python_executable(),
                "-m snakemake",
                format_cli_arg("--snakefile", self.get_snakefile()),
                self.get_job_args(job),
                self.general_args,
                suffix,
            ]
        )


class TouchExecutor(RealExecutor):
    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        try:
            # Touching of output files will be done by handle_job_success
            time.sleep(0.1)
            callback(job)
        except OSError as ex:
            print_exception(ex, self.workflow.linemaps)
            error_callback(job)

    def handle_job_success(self, job):
        super().handle_job_success(job, ignore_missing_output=True)


_ProcessPoolExceptions = (KeyboardInterrupt,)
try:
    from concurrent.futures.process import BrokenProcessPool

    _ProcessPoolExceptions = (KeyboardInterrupt, BrokenProcessPool)
except ImportError:
    pass


class CPUExecutor(RealExecutor):
    def __init__(
        self,
        workflow,
        dag,
        workers,
        printreason=False,
        quiet=False,
        printshellcmds=False,
        use_threads=False,
        cores=1,
        keepincomplete=False,
    ):
        super().__init__(
            workflow,
            dag,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            keepincomplete=keepincomplete,
        )

        self.use_threads = use_threads
        self.cores = cores

        # Zero thread jobs do not need a thread, but they occupy additional workers.
        # Hence we need to reserve additional workers for them.
        self.workers = workers + 5
        self.pool = concurrent.futures.ThreadPoolExecutor(max_workers=self.workers)

    @property
    def job_specific_local_groupid(self):
        return False

    def get_job_exec_prefix(self, job):
        return f"cd {self.workflow.workdir_init}"

    def get_exec_mode(self):
        return Mode.subprocess

    def get_python_executable(self):
        return sys.executable

    def get_envvar_declarations(self):
        return ""

    def get_job_args(self, job, **kwargs):
        return f"{super().get_job_args(job, **kwargs)} --quiet"

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)

        if job.is_group():
            # if we still don't have enough workers for this group, create a new pool here
            missing_workers = max(len(job) - self.workers, 0)
            if missing_workers:
                self.workers += missing_workers
                self.pool = concurrent.futures.ThreadPoolExecutor(
                    max_workers=self.workers
                )

            # the future waits for the entire group job
            future = self.pool.submit(self.run_group_job, job)
        else:
            future = self.run_single_job(job)

        future.add_done_callback(partial(self._callback, job, callback, error_callback))

    def job_args_and_prepare(self, job):
        job.prepare()

        conda_env = (
            job.conda_env.address if self.workflow.use_conda and job.conda_env else None
        )
        container_img = (
            job.container_img_path if self.workflow.use_singularity else None
        )
        env_modules = job.env_modules if self.workflow.use_env_modules else None

        benchmark = None
        benchmark_repeats = job.benchmark_repeats or 1
        if job.benchmark is not None:
            benchmark = str(job.benchmark)
        return (
            job.rule,
            job.input._plainstrings(),
            job.output._plainstrings(),
            job.params,
            job.wildcards,
            job.threads,
            job.resources,
            job.log._plainstrings(),
            benchmark,
            benchmark_repeats,
            conda_env,
            container_img,
            self.workflow.singularity_args,
            env_modules,
            self.workflow.use_singularity,
            self.workflow.linemaps,
            self.workflow.debug,
            self.workflow.cleanup_scripts,
            job.shadow_dir,
            job.jobid,
            self.workflow.edit_notebook if self.dag.is_edit_notebook_job(job) else None,
            self.workflow.conda_base_path,
            job.rule.basedir,
            self.workflow.sourcecache.runtime_cache_path,
        )

    def run_single_job(self, job):
        if (
            self.use_threads
            or (not job.is_shadow and not job.is_run)
            or job.is_template_engine
        ):
            future = self.pool.submit(
                self.cached_or_run, job, run_wrapper, *self.job_args_and_prepare(job)
            )
        else:
            # run directive jobs are spawned into subprocesses
            future = self.pool.submit(self.cached_or_run, job, self.spawn_job, job)
        return future

    def run_group_job(self, job):
        """Run a pipe or service group job.

        This lets all items run simultaneously."""
        # we only have to consider pipe or service groups because in local running mode,
        # these are the only groups that will occur

        futures = [self.run_single_job(j) for j in job]
        n_non_service = sum(1 for j in job if not j.is_service)

        while True:
            n_finished = 0
            for f in futures:
                if f.done():
                    ex = f.exception()
                    if ex is not None:
                        # kill all shell commands of the other group jobs
                        # there can be only shell commands because the
                        # run directive is not allowed for pipe jobs
                        for j in job:
                            shell.kill(j.jobid)
                        raise ex
                    else:
                        n_finished += 1
            if n_finished >= n_non_service:
                # terminate all service jobs since all consumers are done
                for j in job:
                    if j.is_service:
                        logger.info(
                            f"Terminating service job {j.jobid} since all consuming jobs are finished."
                        )
                        shell.terminate(j.jobid)
                        logger.info(
                            f"Service job {j.jobid} has been successfully terminated."
                        )

                return
            time.sleep(1)

    def spawn_job(self, job):
        cmd = self.format_job_exec(job)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            raise SpawnedJobError()

    def cached_or_run(self, job, run_func, *args):
        """
        Either retrieve result from cache, or run job with given function.
        """
        to_cache = self.workflow.is_cached_rule(job.rule)
        try:
            if to_cache:
                self.workflow.output_file_cache.fetch(job)
                return
        except CacheMissException:
            pass
        run_func(*args)
        if to_cache:
            self.workflow.output_file_cache.store(job)

    def shutdown(self):
        self.pool.shutdown()

    def cancel(self):
        self.pool.shutdown()

    def _callback(self, job, callback, error_callback, future):
        try:
            ex = future.exception()
            if ex is not None:
                raise ex
            callback(job)
        except _ProcessPoolExceptions:
            self.handle_job_error(job)
            # no error callback, just silently ignore the interrupt as the main scheduler is also killed
        except SpawnedJobError:
            # don't print error message, this is done by the spawned subprocess
            error_callback(job)
        except (Exception, BaseException) as ex:
            self.print_job_error(job)
            if not (job.is_group() or job.shellcmd) or self.workflow.verbose:
                print_exception(ex, self.workflow.linemaps)
            error_callback(job)

    def handle_job_success(self, job):
        super().handle_job_success(job)

    def handle_job_error(self, job):
        super().handle_job_error(job)
        if not self.keepincomplete:
            job.cleanup()
            self.workflow.persistence.cleanup(job)


class ClusterExecutor(RealExecutor):
    """Backend for distributed execution.

    The key idea is that a job is converted into a script that invokes Snakemake again, in whatever environment is targeted. The script is submitted to some job management platform (e.g. a cluster scheduler like slurm).
    This class can be specialized to generate more specific backends, also for the cloud.
    """

    default_jobscript = "jobscript.sh"

    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob.{name}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        cluster_config=None,
        local_input=None,
        restart_times=None,
        exec_job=None,
        assume_shared_fs=True,
        max_status_checks_per_second=1,
        disable_default_remote_provider_args=False,
        disable_default_resources_args=False,
        disable_envvar_declarations=False,
        keepincomplete=False,
    ):
        from ratelimiter import RateLimiter

        local_input = local_input or []
        super().__init__(
            workflow,
            dag,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            assume_shared_fs=assume_shared_fs,
            keepincomplete=keepincomplete,
        )

        if not self.assume_shared_fs:
            # use relative path to Snakefile
            self.snakefile = os.path.relpath(workflow.main_snakefile)

        self.is_default_jobscript = False
        jobscript = workflow.jobscript
        if jobscript is None:
            jobscript = os.path.join(os.path.dirname(__file__), self.default_jobscript)
            self.is_default_jobscript = True
        try:
            with open(jobscript) as f:
                self.jobscript = f.read()
        except IOError as e:
            raise WorkflowError(e)

        if not "jobid" in get_wildcard_names(jobname):
            raise WorkflowError(
                'Defined jobname ("{}") has to contain the wildcard {jobid}.'
            )

        self.jobname = jobname
        self._tmpdir = None
        self.cores = cores if cores else "all"
        self.cluster_config = cluster_config if cluster_config else dict()

        self.restart_times = restart_times

        self.active_jobs = list()
        self.lock = threading.Lock()
        self.wait = True
        self.wait_thread = threading.Thread(target=self._wait_thread)
        self.wait_thread.daemon = True
        self.wait_thread.start()

        self.disable_default_remote_provider_args = disable_default_remote_provider_args
        self.disable_default_resources_args = disable_default_resources_args
        self.disable_envvar_declarations = disable_envvar_declarations

        self.max_status_checks_per_second = max_status_checks_per_second

        self.status_rate_limiter = RateLimiter(
            max_calls=self.max_status_checks_per_second, period=1
        )

    def get_default_remote_provider_args(self):
        if not self.disable_default_remote_provider_args:
            return super().get_default_remote_provider_args()
        else:
            return ""

    def get_default_resources_args(self, default_resources=None):
        if not self.disable_default_resources_args:
            return super().get_default_resources_args(default_resources)
        else:
            return ""

    def get_workdir_arg(self):
        if self.assume_shared_fs:
            return super().get_workdir_arg()
        return ""

    def get_envvar_declarations(self):
        if not self.disable_envvar_declarations:
            return " ".join(
                f"{var}={repr(os.environ[var])}" for var in self.workflow.envvars
            )
        else:
            return ""

    def get_python_executable(self):
        return sys.executable if self.assume_shared_fs else "python"

    def get_exec_mode(self):
        return Mode.cluster

    def get_job_args(self, job):
        waitfiles_parameter = ""
        if self.assume_shared_fs:
            wait_for_files = []
            wait_for_files.append(self.tmpdir)
            wait_for_files.extend(job.get_wait_for_files())

            # Only create extra file if we have more than 20 input files.
            # This should not require the file creation in most cases.
            if len(wait_for_files) > 20:
                wait_for_files_file = self.get_jobscript(job) + ".waitforfilesfile.txt"
                with open(wait_for_files_file, "w") as fd:
                    print(*wait_for_files, sep="\n", file=fd)

                waitfiles_parameter = format_cli_arg(
                    "--wait-for-files-file", wait_for_files_file
                )
            else:
                waitfiles_parameter = format_cli_arg("--wait-for-files", wait_for_files)

        return f"{super().get_job_args(job)} {waitfiles_parameter}"

    def _wait_thread(self):
        try:
            self._wait_for_jobs()
        except Exception as e:
            self.workflow.scheduler.executor_error_callback(e)

    def shutdown(self):
        with self.lock:
            self.wait = False
        self.wait_thread.join()
        if not self.workflow.immediate_submit:
            # Only delete tmpdir (containing jobscripts) if not using
            # immediate_submit. With immediate_submit, jobs can be scheduled
            # after this method is completed. Hence we have to keep the
            # directory.
            shutil.rmtree(self.tmpdir)

    def cancel(self):
        self.shutdown()

    def _run(self, job, callback=None, error_callback=None):
        if self.assume_shared_fs:
            job.remove_existing_output()
            job.download_remote_input()
        super()._run(job, callback=callback, error_callback=error_callback)

    @property
    def tmpdir(self):
        if self._tmpdir is None:
            self._tmpdir = tempfile.mkdtemp(dir=".snakemake", prefix="tmp.")
        return os.path.abspath(self._tmpdir)

    def get_jobscript(self, job):
        f = job.format_wildcards(self.jobname, cluster=self.cluster_wildcards(job))

        if os.path.sep in f:
            raise WorkflowError(
                "Path separator ({}) found in job name {}. "
                "This is not supported.".format(os.path.sep, f)
            )

        return os.path.join(self.tmpdir, f)

    def write_jobscript(self, job, jobscript):
        exec_job = self.format_job_exec(job)

        try:
            content = self.jobscript.format(
                properties=job.properties(cluster=self.cluster_params(job)),
                exec_job=exec_job,
            )
        except KeyError as e:
            if self.is_default_jobscript:
                raise e
            else:
                raise WorkflowError(
                    f"Error formatting custom jobscript {self.workflow.jobscript}: value for {e} not found.\n"
                    "Make sure that your custom jobscript is defined as expected."
                )

        logger.debug("Jobscript:\n{}".format(content))
        with open(jobscript, "w") as f:
            print(content, file=f)
        os.chmod(jobscript, os.stat(jobscript).st_mode | stat.S_IXUSR | stat.S_IRUSR)

    def cluster_params(self, job):
        """Return wildcards object for job from cluster_config."""
        cluster = self.cluster_config.get("__default__", dict()).copy()
        cluster.update(self.cluster_config.get(job.name, dict()))
        # Format values with available parameters from the job.
        for key, value in list(cluster.items()):
            if isinstance(value, str):
                try:
                    cluster[key] = job.format_wildcards(value)
                except NameError as e:
                    if job.is_group():
                        msg = (
                            "Failed to format cluster config for group job. "
                            "You have to ensure that your default entry "
                            "does not contain any items that group jobs "
                            "cannot provide, like {rule}, {wildcards}."
                        )
                    else:
                        msg = (
                            "Failed to format cluster config "
                            "entry for job {}.".format(job.rule.name)
                        )
                    raise WorkflowError(msg, e)

        return cluster

    def cluster_wildcards(self, job):
        return Wildcards(fromdict=self.cluster_params(job))

    def handle_job_success(self, job):
        super().handle_job_success(
            job, upload_remote=False, handle_log=False, handle_touch=False
        )

    def handle_job_error(self, job):
        # TODO what about removing empty remote dirs?? This cannot be decided
        # on the cluster node.
        super().handle_job_error(job, upload_remote=False)
        logger.debug("Cleanup job metadata.")
        # We have to remove metadata here as well.
        # It will be removed by the CPUExecutor in case of a shared FS,
        # but we might not see the removal due to filesystem latency.
        # By removing it again, we make sure that it is gone on the host FS.
        if not self.keepincomplete:
            self.workflow.persistence.cleanup(job)
            # Also cleanup the jobs output files, in case the remote job
            # was not able to, due to e.g. timeout.
            logger.debug("Cleanup failed jobs output files.")
            job.cleanup()

    def print_cluster_job_error(self, job_info, jobid):
        job = job_info.job
        kind = (
            "rule {}".format(job.rule.name)
            if not job.is_group()
            else "group job {}".format(job.groupid)
        )
        logger.error(
            "Error executing {} on cluster (jobid: {}, external: "
            "{}, jobscript: {}). For error details see the cluster "
            "log and the log files of the involved rule(s).".format(
                kind, jobid, job_info.jobid, job_info.jobscript
            )
        )


GenericClusterJob = namedtuple(
    "GenericClusterJob",
    "job jobid callback error_callback jobscript jobfinished jobfailed",
)


class GenericClusterExecutor(ClusterExecutor):
    def __init__(
        self,
        workflow,
        dag,
        cores,
        submitcmd="qsub",
        statuscmd=None,
        cancelcmd=None,
        cancelnargs=None,
        sidecarcmd=None,
        cluster_config=None,
        jobname="snakejob.{rulename}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        restart_times=0,
        assume_shared_fs=True,
        max_status_checks_per_second=1,
        keepincomplete=False,
    ):

        self.submitcmd = submitcmd
        if not assume_shared_fs and statuscmd is None:
            raise WorkflowError(
                "When no shared filesystem can be assumed, a "
                "status command must be given."
            )

        self.statuscmd = statuscmd
        self.cancelcmd = cancelcmd
        self.sidecarcmd = sidecarcmd
        self.cancelnargs = cancelnargs
        self.external_jobid = dict()
        # We need to collect all external ids so we can properly cancel even if
        # the status update queue is running.
        self.all_ext_jobids = list()

        super().__init__(
            workflow,
            dag,
            cores,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            cluster_config=cluster_config,
            restart_times=restart_times,
            assume_shared_fs=assume_shared_fs,
            max_status_checks_per_second=max_status_checks_per_second,
            keepincomplete=keepincomplete,
        )

        self.sidecar_vars = None
        if self.sidecarcmd:
            self._launch_sidecar()

        if not statuscmd and not assume_shared_fs:
            raise WorkflowError(
                "If no shared filesystem is used, you have to "
                "specify a cluster status command."
            )

    def get_job_exec_prefix(self, job):
        if self.assume_shared_fs:
            # quoting the workdir since it may contain spaces
            return f"cd {repr(self.workflow.workdir_init)}"
        else:
            return ""

    def get_job_exec_suffix(self, job):
        if self.statuscmd:
            return "exit 0 || exit 1"
        elif self.assume_shared_fs:
            # TODO wrap with watch and touch {jobrunning}
            # check modification date of {jobrunning} in the wait_for_job method

            return (
                f"touch {repr(self.get_jobfinished_marker(job))} || "
                f"(touch {repr(self.get_jobfailed_marker(job))}; exit 1)"
            )
        assert False, "bug: neither statuscmd defined nor shared FS"

    def get_jobfinished_marker(self, job):
        return os.path.join(self.tmpdir, "{}.jobfinished".format(job.jobid))

    def get_jobfailed_marker(self, job):
        return os.path.join(self.tmpdir, "{}.jobfailed".format(job.jobid))

    def _launch_sidecar(self):
        def copy_stdout(executor, process):
            """Run sidecar process and copy it's stdout to our stdout."""
            while process.poll() is None and executor.wait:
                buf = process.stdout.readline()
                if buf:
                    sys.stdout.write(buf)
            # one final time ...
            buf = process.stdout.readline()
            if buf:
                sys.stdout.write(buf)

        def wait(executor, process):
            while executor.wait:
                time.sleep(0.5)
            process.terminate()
            process.wait()
            logger.info(
                "Cluster sidecar process has terminated (retcode=%d)."
                % process.returncode
            )

        logger.info("Launch sidecar process and read first output line.")
        process = subprocess.Popen(
            self.sidecarcmd, stdout=subprocess.PIPE, shell=False, encoding="utf-8"
        )
        self.sidecar_vars = process.stdout.readline()
        while self.sidecar_vars and self.sidecar_vars[-1] in "\n\r":
            self.sidecar_vars = self.sidecar_vars[:-1]
        logger.info("Done reading first output line.")

        thread_stdout = threading.Thread(
            target=copy_stdout, name="sidecar_stdout", args=(self, process)
        )
        thread_stdout.start()
        thread_wait = threading.Thread(
            target=wait, name="sidecar_stdout", args=(self, process)
        )
        thread_wait.start()

    def cancel(self):
        def _chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            for i in range(0, len(lst), n):
                yield lst[i : i + n]

        if self.cancelcmd:  # We have --cluster-cancel
            # Enumerate job IDs and create chunks.  If cancelnargs evaluates to false (0/None)
            # then pass all job ids at once
            jobids = list(self.all_ext_jobids)
            chunks = list(_chunks(jobids, self.cancelnargs or len(jobids)))
            # Go through the chunks and cancel the jobs, warn in case of failures.
            failures = 0
            for chunk in chunks:
                try:
                    cancel_timeout = 2  # rather fail on timeout than miss canceling all
                    env = dict(os.environ)
                    if self.sidecar_vars:
                        env["SNAKEMAKE_CLUSTER_SIDECAR_VARS"] = self.sidecar_vars
                    subprocess.check_call(
                        [self.cancelcmd] + chunk,
                        shell=False,
                        timeout=cancel_timeout,
                        env=env,
                    )
                except subprocess.SubprocessError:
                    failures += 1
            if failures:
                logger.info(
                    (
                        "{} out of {} calls to --cluster-cancel failed.  This is safe to "
                        "ignore in most cases."
                    ).format(failures, len(chunks))
                )
        else:
            logger.info(
                "No --cluster-cancel given. Will exit after finishing currently running jobs."
            )
            self.shutdown()

    def register_job(self, job):
        # Do not register job here.
        # Instead do it manually once the jobid is known.
        pass

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        jobid = job.jobid

        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        jobfinished = self.get_jobfinished_marker(job)
        jobfailed = self.get_jobfailed_marker(job)

        if self.statuscmd:
            ext_jobid = self.dag.incomplete_external_jobid(job)
            if ext_jobid:
                # Job is incomplete and still running.
                # We simply register it and wait for completion or failure.
                logger.info(
                    "Resuming incomplete job {} with external jobid '{}'.".format(
                        jobid, ext_jobid
                    )
                )
                submit_callback(job)
                with self.lock:
                    self.all_ext_jobids.append(ext_jobid)
                    self.active_jobs.append(
                        GenericClusterJob(
                            job,
                            ext_jobid,
                            callback,
                            error_callback,
                            jobscript,
                            jobfinished,
                            jobfailed,
                        )
                    )
                return

        deps = " ".join(
            self.external_jobid[f] for f in job.input if f in self.external_jobid
        )
        try:
            submitcmd = job.format_wildcards(
                self.submitcmd, dependencies=deps, cluster=self.cluster_wildcards(job)
            )
        except AttributeError as e:
            raise WorkflowError(str(e), rule=job.rule if not job.is_group() else None)

        try:
            env = dict(os.environ)
            if self.sidecar_vars:
                env["SNAKEMAKE_CLUSTER_SIDECAR_VARS"] = self.sidecar_vars

            # Remove SNAKEMAKE_PROFILE from environment as the snakemake call inside
            # of the cluster job must run locally (or complains about missing -j).
            env.pop("SNAKEMAKE_PROFILE", None)

            ext_jobid = (
                subprocess.check_output(
                    '{submitcmd} "{jobscript}"'.format(
                        submitcmd=submitcmd, jobscript=jobscript
                    ),
                    shell=True,
                    env=env,
                )
                .decode()
                .split("\n")
            )
        except subprocess.CalledProcessError as ex:
            logger.error(
                "Error submitting jobscript (exit code {}):\n{}".format(
                    ex.returncode, ex.output.decode()
                )
            )
            error_callback(job)
            return
        if ext_jobid and ext_jobid[0]:
            ext_jobid = ext_jobid[0]
            self.external_jobid.update((f, ext_jobid) for f in job.output)
            logger.info(
                "Submitted {} {} with external jobid '{}'.".format(
                    "group job" if job.is_group() else "job", jobid, ext_jobid
                )
            )
            self.workflow.persistence.started(job, external_jobid=ext_jobid)

        submit_callback(job)

        with self.lock:
            self.all_ext_jobids.append(ext_jobid)
            self.active_jobs.append(
                GenericClusterJob(
                    job,
                    ext_jobid,
                    callback,
                    error_callback,
                    jobscript,
                    jobfinished,
                    jobfailed,
                )
            )

    def _wait_for_jobs(self):
        success = "success"
        failed = "failed"
        running = "running"
        status_cmd_kills = []
        if self.statuscmd is not None:

            def job_status(job, valid_returns=["running", "success", "failed"]):
                try:
                    # this command shall return "success", "failed" or "running"
                    env = dict(os.environ)
                    if self.sidecar_vars:
                        env["SNAKEMAKE_CLUSTER_SIDECAR_VARS"] = self.sidecar_vars
                    ret = subprocess.check_output(
                        "{statuscmd} '{jobid}'".format(
                            jobid=job.jobid, statuscmd=self.statuscmd
                        ),
                        shell=True,
                        env=env,
                    ).decode()
                except subprocess.CalledProcessError as e:
                    if e.returncode < 0:
                        # Ignore SIGINT and all other issues due to signals
                        # because it will be caused by hitting e.g.
                        # Ctrl-C on the main process or sending killall to
                        # snakemake.
                        # Snakemake will handle the signal in
                        # the main process.
                        status_cmd_kills.append(-e.returncode)
                        if len(status_cmd_kills) > 10:
                            logger.info(
                                "Cluster status command {} was killed >10 times with signal(s) {} "
                                "(if this happens unexpectedly during your workflow execution, "
                                "have a closer look.).".format(
                                    self.statuscmd, ",".join(status_cmd_kills)
                                )
                            )
                            status_cmd_kills.clear()
                    else:
                        raise WorkflowError(
                            "Failed to obtain job status. "
                            "See above for error message."
                        )

                ret = ret.strip().split("\n")
                if len(ret) != 1 or ret[0] not in valid_returns:
                    raise WorkflowError(
                        "Cluster status command {} returned {} but just a single line with one of {} is expected.".format(
                            self.statuscmd, "\\n".join(ret), ",".join(valid_returns)
                        )
                    )
                return ret[0]

        else:

            def job_status(job):
                if os.path.exists(active_job.jobfinished):
                    os.remove(active_job.jobfinished)
                    os.remove(active_job.jobscript)
                    return success
                if os.path.exists(active_job.jobfailed):
                    os.remove(active_job.jobfailed)
                    os.remove(active_job.jobscript)
                    return failed
                return running

        while True:
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            # logger.debug("Checking status of {} jobs.".format(len(active_jobs)))
            for active_job in active_jobs:
                with self.status_rate_limiter:
                    status = job_status(active_job)

                    if status == success:
                        active_job.callback(active_job.job)
                    elif status == failed:
                        self.print_job_error(
                            active_job.job,
                            cluster_jobid=active_job.jobid
                            if active_job.jobid
                            else "unknown",
                        )
                        self.print_cluster_job_error(
                            active_job, self.dag.jobid(active_job.job)
                        )
                        active_job.error_callback(active_job.job)
                    else:
                        still_running.append(active_job)
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


SynchronousClusterJob = namedtuple(
    "SynchronousClusterJob", "job jobid callback error_callback jobscript process"
)


class SynchronousClusterExecutor(ClusterExecutor):
    """
    invocations like "qsub -sync y" (SGE) or "bsub -K" (LSF) are
    synchronous, blocking the foreground thread and returning the
    remote exit code at remote exit.
    """

    def __init__(
        self,
        workflow,
        dag,
        cores,
        submitcmd="qsub",
        cluster_config=None,
        jobname="snakejob.{rulename}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        restart_times=0,
        assume_shared_fs=True,
        keepincomplete=False,
    ):
        super().__init__(
            workflow,
            dag,
            cores,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            cluster_config=cluster_config,
            restart_times=restart_times,
            assume_shared_fs=assume_shared_fs,
            max_status_checks_per_second=10,
            keepincomplete=keepincomplete,
        )
        self.submitcmd = submitcmd
        self.external_jobid = dict()

    def get_job_exec_prefix(self, job):
        if self.assume_shared_fs:
            return f"cd {self.workflow.workdir_init}"
        else:
            return ""

    def cancel(self):
        logger.info("Will exit after finishing currently running jobs.")
        self.shutdown()

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)

        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        deps = " ".join(
            self.external_jobid[f] for f in job.input if f in self.external_jobid
        )
        try:
            submitcmd = job.format_wildcards(
                self.submitcmd, dependencies=deps, cluster=self.cluster_wildcards(job)
            )
        except AttributeError as e:
            raise WorkflowError(str(e), rule=job.rule if not job.is_group() else None)

        process = subprocess.Popen(
            '{submitcmd} "{jobscript}"'.format(
                submitcmd=submitcmd, jobscript=jobscript
            ),
            shell=True,
        )
        submit_callback(job)

        with self.lock:
            self.active_jobs.append(
                SynchronousClusterJob(
                    job, process.pid, callback, error_callback, jobscript, process
                )
            )

    def _wait_for_jobs(self):
        while True:
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for active_job in active_jobs:
                with self.status_rate_limiter:
                    exitcode = active_job.process.poll()
                    if exitcode is None:
                        # job not yet finished
                        still_running.append(active_job)
                    elif exitcode == 0:
                        # job finished successfully
                        os.remove(active_job.jobscript)
                        active_job.callback(active_job.job)
                    else:
                        # job failed
                        os.remove(active_job.jobscript)
                        self.print_job_error(active_job.job)
                        self.print_cluster_job_error(
                            active_job, self.dag.jobid(active_job.job)
                        )
                        active_job.error_callback(active_job.job)
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


DRMAAClusterJob = namedtuple(
    "DRMAAClusterJob", "job jobid callback error_callback jobscript"
)


class DRMAAExecutor(ClusterExecutor):
    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob.{rulename}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        drmaa_args="",
        drmaa_log_dir=None,
        cluster_config=None,
        restart_times=0,
        assume_shared_fs=True,
        max_status_checks_per_second=1,
        keepincomplete=False,
    ):
        super().__init__(
            workflow,
            dag,
            cores,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            cluster_config=cluster_config,
            restart_times=restart_times,
            assume_shared_fs=assume_shared_fs,
            max_status_checks_per_second=max_status_checks_per_second,
            keepincomplete=keepincomplete,
        )
        try:
            import drmaa
        except ImportError:
            raise WorkflowError(
                "Python support for DRMAA is not installed. "
                "Please install it, e.g. with easy_install3 --user drmaa"
            )
        except RuntimeError as e:
            raise WorkflowError("Error loading drmaa support:\n{}".format(e))
        self.session = drmaa.Session()
        self.drmaa_args = drmaa_args
        self.drmaa_log_dir = drmaa_log_dir
        self.session.initialize()
        self.submitted = list()

    def get_job_exec_prefix(self, job):
        if self.assume_shared_fs:
            return f"cd {self.workflow.workdir_init}"
        else:
            return ""

    def cancel(self):
        from drmaa.const import JobControlAction
        from drmaa.errors import InvalidJobException, InternalException

        for jobid in self.submitted:
            try:
                self.session.control(jobid, JobControlAction.TERMINATE)
            except (InvalidJobException, InternalException):
                # This is common - logging a warning would probably confuse the user.
                pass
        self.shutdown()

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        try:
            drmaa_args = job.format_wildcards(
                self.drmaa_args, cluster=self.cluster_wildcards(job)
            )
        except AttributeError as e:
            raise WorkflowError(str(e), rule=job.rule)

        import drmaa

        if self.drmaa_log_dir:
            makedirs(self.drmaa_log_dir)

        try:
            jt = self.session.createJobTemplate()
            jt.remoteCommand = jobscript
            jt.nativeSpecification = drmaa_args
            if self.drmaa_log_dir:
                jt.outputPath = ":" + self.drmaa_log_dir
                jt.errorPath = ":" + self.drmaa_log_dir
            jt.jobName = os.path.basename(jobscript)

            jobid = self.session.runJob(jt)
        except (
            drmaa.DeniedByDrmException,
            drmaa.InternalException,
            drmaa.InvalidAttributeValueException,
        ) as e:
            print_exception(
                WorkflowError("DRMAA Error: {}".format(e)), self.workflow.linemaps
            )
            error_callback(job)
            return
        logger.info(
            "Submitted DRMAA job {} with external jobid {}.".format(job.jobid, jobid)
        )
        self.submitted.append(jobid)
        self.session.deleteJobTemplate(jt)

        submit_callback(job)

        with self.lock:
            self.active_jobs.append(
                DRMAAClusterJob(job, jobid, callback, error_callback, jobscript)
            )

    def shutdown(self):
        super().shutdown()
        self.session.exit()

    def _wait_for_jobs(self):
        import drmaa

        suspended_msg = set()

        while True:
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for active_job in active_jobs:
                with self.status_rate_limiter:
                    try:
                        retval = self.session.jobStatus(active_job.jobid)
                    except drmaa.ExitTimeoutException as e:
                        # job still active
                        still_running.append(active_job)
                        continue
                    except (drmaa.InternalException, Exception) as e:
                        print_exception(
                            WorkflowError("DRMAA Error: {}".format(e)),
                            self.workflow.linemaps,
                        )
                        os.remove(active_job.jobscript)
                        active_job.error_callback(active_job.job)
                        continue
                    if retval == drmaa.JobState.DONE:
                        os.remove(active_job.jobscript)
                        active_job.callback(active_job.job)
                    elif retval == drmaa.JobState.FAILED:
                        os.remove(active_job.jobscript)
                        self.print_job_error(active_job.job)
                        self.print_cluster_job_error(
                            active_job, self.dag.jobid(active_job.job)
                        )
                        active_job.error_callback(active_job.job)
                    else:
                        # still running
                        still_running.append(active_job)

                        def handle_suspended(by):
                            if active_job.job.jobid not in suspended_msg:
                                logger.warning(
                                    "Job {} (DRMAA id: {}) was suspended by {}.".format(
                                        active_job.job.jobid, active_job.jobid, by
                                    )
                                )
                                suspended_msg.add(active_job.job.jobid)

                        if retval == drmaa.JobState.USER_SUSPENDED:
                            handle_suspended("user")
                        elif retval == drmaa.JobState.SYSTEM_SUSPENDED:
                            handle_suspended("system")
                        else:
                            try:
                                suspended_msg.remove(active_job.job.jobid)
                            except KeyError:
                                # there was nothing to remove
                                pass

            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


@contextlib.contextmanager
def change_working_directory(directory=None):
    """Change working directory in execution context if provided."""
    if directory:
        try:
            saved_directory = os.getcwd()
            logger.info("Changing to shadow directory: {}".format(directory))
            os.chdir(directory)
            yield
        finally:
            os.chdir(saved_directory)
    else:
        yield


KubernetesJob = namedtuple(
    "KubernetesJob", "job jobid callback error_callback kubejob jobscript"
)


class KubernetesExecutor(ClusterExecutor):
    def __init__(
        self,
        workflow,
        dag,
        namespace,
        container_image=None,
        jobname="{rulename}.{jobid}",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        cluster_config=None,
        local_input=None,
        restart_times=None,
        keepincomplete=False,
    ):
        self.workflow = workflow

        super().__init__(
            workflow,
            dag,
            None,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            cluster_config=cluster_config,
            local_input=local_input,
            restart_times=restart_times,
            assume_shared_fs=False,
            max_status_checks_per_second=10,
            disable_envvar_declarations=True,
        )
        # use relative path to Snakefile
        self.snakefile = os.path.relpath(workflow.main_snakefile)

        try:
            from kubernetes import config
        except ImportError:
            raise WorkflowError(
                "The Python 3 package 'kubernetes' "
                "must be installed to use Kubernetes"
            )
        config.load_kube_config()

        import kubernetes.client

        self.kubeapi = kubernetes.client.CoreV1Api()
        self.batchapi = kubernetes.client.BatchV1Api()
        self.namespace = namespace
        self.envvars = workflow.envvars
        self.secret_files = {}
        self.run_namespace = str(uuid.uuid4())
        self.secret_envvars = {}
        self.register_secret()
        self.container_image = container_image or get_container_image()
        logger.info(f"Using {self.container_image} for Kubernetes jobs.")

    def get_job_exec_prefix(self, job):
        return "cp -rf /source/. ."

    def register_secret(self):
        import kubernetes.client

        secret = kubernetes.client.V1Secret()
        secret.metadata = kubernetes.client.V1ObjectMeta()
        # create a random uuid
        secret.metadata.name = self.run_namespace
        secret.type = "Opaque"
        secret.data = {}
        for i, f in enumerate(self.workflow.get_sources()):
            if f.startswith(".."):
                logger.warning(
                    "Ignoring source file {}. Only files relative "
                    "to the working directory are allowed.".format(f)
                )
                continue

            # The kubernetes API can't create secret files larger than 1MB.
            source_file_size = os.path.getsize(f)
            max_file_size = 1048576
            if source_file_size > max_file_size:
                logger.warning(
                    "Skipping the source file {f}. Its size {source_file_size} exceeds "
                    "the maximum file size (1MB) that can be passed "
                    "from host to kubernetes.".format(
                        f=f, source_file_size=source_file_size
                    )
                )
                continue

            with open(f, "br") as content:
                key = "f{}".format(i)

                # Some files are smaller than 1MB, but grows larger after being base64 encoded
                # We should exclude them as well, otherwise Kubernetes APIs will complain
                encoded_contents = base64.b64encode(content.read()).decode()
                encoded_size = len(encoded_contents)
                if encoded_size > 1048576:
                    logger.warning(
                        "Skipping the source file {f} for secret key {key}. "
                        "Its base64 encoded size {encoded_size} exceeds "
                        "the maximum file size (1MB) that can be passed "
                        "from host to kubernetes.".format(
                            f=f,
                            source_file_size=source_file_size,
                            key=key,
                            encoded_size=encoded_size,
                        )
                    )
                    continue

                self.secret_files[key] = f
                secret.data[key] = encoded_contents

        for e in self.envvars:
            try:
                key = e.lower()
                secret.data[key] = base64.b64encode(os.environ[e].encode()).decode()
                self.secret_envvars[key] = e
            except KeyError:
                continue

        # Test if the total size of the configMap exceeds 1MB
        config_map_size = sum(
            [len(base64.b64decode(v)) for k, v in secret.data.items()]
        )
        if config_map_size > 1048576:
            logger.warning(
                "The total size of the included files and other Kubernetes secrets "
                "is {}, exceeding the 1MB limit.\n".format(config_map_size)
            )
            logger.warning(
                "The following are the largest files. Consider removing some of them "
                "(you need remove at least {} bytes):".format(config_map_size - 1048576)
            )

            entry_sizes = {
                self.secret_files[k]: len(base64.b64decode(v))
                for k, v in secret.data.items()
                if k in self.secret_files
            }
            for k, v in sorted(entry_sizes.items(), key=lambda item: item[1])[:-6:-1]:
                logger.warning("  * File: {k}, original size: {v}".format(k=k, v=v))

            raise WorkflowError("ConfigMap too large")

        self.kubeapi.create_namespaced_secret(self.namespace, secret)

    def unregister_secret(self):
        import kubernetes.client

        safe_delete_secret = lambda: self.kubeapi.delete_namespaced_secret(
            self.run_namespace, self.namespace, body=kubernetes.client.V1DeleteOptions()
        )
        self._kubernetes_retry(safe_delete_secret)

    # In rare cases, deleting a pod may rais 404 NotFound error.
    def safe_delete_pod(self, jobid, ignore_not_found=True):
        import kubernetes.client

        body = kubernetes.client.V1DeleteOptions()
        try:
            self.kubeapi.delete_namespaced_pod(jobid, self.namespace, body=body)
        except kubernetes.client.rest.ApiException as e:
            if e.status == 404 and ignore_not_found:
                # Can't find the pod. Maybe it's already been
                # destroyed. Proceed with a warning message.
                logger.warning(
                    "[WARNING] 404 not found when trying to delete the pod: {jobid}\n"
                    "[WARNING] Ignore this error\n".format(jobid=jobid)
                )
            else:
                raise e

    def shutdown(self):
        self.unregister_secret()
        super().shutdown()

    def cancel(self):
        import kubernetes.client

        body = kubernetes.client.V1DeleteOptions()
        with self.lock:
            for j in self.active_jobs:
                func = lambda: self.safe_delete_pod(j.jobid, ignore_not_found=True)
                self._kubernetes_retry(func)

        self.shutdown()

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        import kubernetes.client

        super()._run(job)
        exec_job = self.format_job_exec(job)

        # Kubernetes silently does not submit a job if the name is too long
        # therefore, we ensure that it is not longer than snakejob+uuid.
        jobid = "snakejob-{}".format(
            get_uuid("{}-{}-{}".format(self.run_namespace, job.jobid, job.attempt))
        )

        body = kubernetes.client.V1Pod()
        body.metadata = kubernetes.client.V1ObjectMeta(labels={"app": "snakemake"})

        body.metadata.name = jobid

        # container
        container = kubernetes.client.V1Container(name=jobid)
        container.image = self.container_image
        container.command = shlex.split("/bin/sh")
        container.args = ["-c", exec_job]
        container.working_dir = "/workdir"
        container.volume_mounts = [
            kubernetes.client.V1VolumeMount(name="workdir", mount_path="/workdir")
        ]
        container.volume_mounts = [
            kubernetes.client.V1VolumeMount(name="source", mount_path="/source")
        ]

        node_selector = {}
        if "machine_type" in job.resources.keys():
            # Kubernetes labels a node by its instance type using this node_label.
            node_selector["node.kubernetes.io/instance-type"] = job.resources[
                "machine_type"
            ]

        body.spec = kubernetes.client.V1PodSpec(
            containers=[container], node_selector=node_selector
        )
        # fail on first error
        body.spec.restart_policy = "Never"

        # source files as a secret volume
        # we copy these files to the workdir before executing Snakemake
        too_large = [
            path
            for path in self.secret_files.values()
            if os.path.getsize(path) > 1000000
        ]
        if too_large:
            raise WorkflowError(
                "The following source files exceed the maximum "
                "file size (1MB) that can be passed from host to "
                "kubernetes. These are likely not source code "
                "files. Consider adding them to your "
                "remote storage instead or (if software) use "
                "Conda packages or container images:\n{}".format("\n".join(too_large))
            )
        secret_volume = kubernetes.client.V1Volume(name="source")
        secret_volume.secret = kubernetes.client.V1SecretVolumeSource()
        secret_volume.secret.secret_name = self.run_namespace
        secret_volume.secret.items = [
            kubernetes.client.V1KeyToPath(key=key, path=path)
            for key, path in self.secret_files.items()
        ]
        # workdir as an emptyDir volume of undefined size
        workdir_volume = kubernetes.client.V1Volume(name="workdir")
        workdir_volume.empty_dir = kubernetes.client.V1EmptyDirVolumeSource()
        body.spec.volumes = [secret_volume, workdir_volume]

        # env vars
        container.env = []
        for key, e in self.secret_envvars.items():
            envvar = kubernetes.client.V1EnvVar(name=e)
            envvar.value_from = kubernetes.client.V1EnvVarSource()
            envvar.value_from.secret_key_ref = kubernetes.client.V1SecretKeySelector(
                key=key, name=self.run_namespace
            )
            container.env.append(envvar)

        # request resources
        container.resources = kubernetes.client.V1ResourceRequirements()
        container.resources.requests = {}
        container.resources.requests["cpu"] = job.resources["_cores"]
        if "mem_mb" in job.resources.keys():
            container.resources.requests["memory"] = "{}M".format(
                job.resources["mem_mb"]
            )

        # capabilities
        if job.needs_singularity and self.workflow.use_singularity:
            # TODO this should work, but it doesn't currently because of
            # missing loop devices
            # singularity inside docker requires SYS_ADMIN capabilities
            # see https://groups.google.com/a/lbl.gov/forum/#!topic/singularity/e9mlDuzKowc
            # container.capabilities = kubernetes.client.V1Capabilities()
            # container.capabilities.add = ["SYS_ADMIN",
            #                               "DAC_OVERRIDE",
            #                               "SETUID",
            #                               "SETGID",
            #                               "SYS_CHROOT"]

            # Running in priviledged mode always works
            container.security_context = kubernetes.client.V1SecurityContext(
                privileged=True
            )

        pod = self._kubernetes_retry(
            lambda: self.kubeapi.create_namespaced_pod(self.namespace, body)
        )

        logger.info(
            "Get status with:\n"
            "kubectl describe pod {jobid}\n"
            "kubectl logs {jobid}".format(jobid=jobid)
        )
        self.active_jobs.append(
            KubernetesJob(job, jobid, callback, error_callback, pod, None)
        )

    # Sometimes, certain k8s requests throw kubernetes.client.rest.ApiException
    # Solving this issue requires reauthentication, as _kubernetes_retry shows
    # However, reauthentication itself, under rare conditions, may also throw
    # errors such as:
    #   kubernetes.client.exceptions.ApiException: (409), Reason: Conflict
    #
    # This error doesn't mean anything wrong with the k8s cluster, and users can safely
    # ignore it.
    def _reauthenticate_and_retry(self, func=None):
        import kubernetes

        # Unauthorized.
        # Reload config in order to ensure token is
        # refreshed. Then try again.
        logger.info("Trying to reauthenticate")
        kubernetes.config.load_kube_config()
        subprocess.run(["kubectl", "get", "nodes"])

        self.kubeapi = kubernetes.client.CoreV1Api()
        self.batchapi = kubernetes.client.BatchV1Api()

        try:
            self.register_secret()
        except kubernetes.client.rest.ApiException as e:
            if e.status == 409 and e.reason == "Conflict":
                logger.warning("409 conflict ApiException when registering secrets")
                logger.warning(e)
            else:
                raise WorkflowError(
                    e,
                    "This is likely a bug in "
                    "https://github.com/kubernetes-client/python.",
                )

        if func:
            return func()

    def _kubernetes_retry(self, func):
        import kubernetes
        import urllib3

        with self.lock:
            try:
                return func()
            except kubernetes.client.rest.ApiException as e:
                if e.status == 401:
                    # Unauthorized.
                    # Reload config in order to ensure token is
                    # refreshed. Then try again.
                    return self._reauthenticate_and_retry(func)
            # Handling timeout that may occur in case of GKE master upgrade
            except urllib3.exceptions.MaxRetryError as e:
                logger.info(
                    "Request time out! "
                    "check your connection to Kubernetes master"
                    "Workflow will pause for 5 minutes to allow any update operations to complete"
                )
                time.sleep(300)
                try:
                    return func()
                except:
                    # Still can't reach the server after 5 minutes
                    raise WorkflowError(
                        e,
                        "Error 111 connection timeout, please check"
                        " that the k8 cluster master is reachable!",
                    )

    def _wait_for_jobs(self):
        import kubernetes

        while True:
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for j in active_jobs:
                with self.status_rate_limiter:
                    logger.debug("Checking status for pod {}".format(j.jobid))
                    job_not_found = False
                    try:
                        res = self._kubernetes_retry(
                            lambda: self.kubeapi.read_namespaced_pod_status(
                                j.jobid, self.namespace
                            )
                        )
                    except kubernetes.client.rest.ApiException as e:
                        if e.status == 404:
                            # Jobid not found
                            # The job is likely already done and was deleted on
                            # the server.
                            j.callback(j.job)
                            continue
                    except WorkflowError as e:
                        print_exception(e, self.workflow.linemaps)
                        j.error_callback(j.job)
                        continue

                    if res is None:
                        msg = (
                            "Unknown pod {jobid}. "
                            "Has the pod been deleted "
                            "manually?"
                        ).format(jobid=j.jobid)
                        self.print_job_error(j.job, msg=msg, jobid=j.jobid)
                        j.error_callback(j.job)
                    elif res.status.phase == "Failed":
                        msg = (
                            "For details, please issue:\n"
                            "kubectl describe pod {jobid}\n"
                            "kubectl logs {jobid}"
                        ).format(jobid=j.jobid)
                        # failed
                        self.print_job_error(j.job, msg=msg, jobid=j.jobid)
                        j.error_callback(j.job)
                    elif res.status.phase == "Succeeded":
                        # finished
                        j.callback(j.job)

                        func = lambda: self.safe_delete_pod(
                            j.jobid, ignore_not_found=True
                        )
                        self._kubernetes_retry(func)
                    else:
                        # still active
                        still_running.append(j)
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


TibannaJob = namedtuple(
    "TibannaJob", "job jobname jobid exec_arn callback error_callback"
)


class TibannaExecutor(ClusterExecutor):
    def __init__(
        self,
        workflow,
        dag,
        cores,
        tibanna_sfn,
        precommand="",
        tibanna_config=False,
        container_image=None,
        printreason=False,
        quiet=False,
        printshellcmds=False,
        local_input=None,
        restart_times=None,
        max_status_checks_per_second=1,
        keepincomplete=False,
    ):
        self.workflow = workflow
        self.workflow_sources = []
        for wfs in workflow.get_sources():
            if os.path.isdir(wfs):
                for (dirpath, dirnames, filenames) in os.walk(wfs):
                    self.workflow_sources.extend(
                        [os.path.join(dirpath, f) for f in filenames]
                    )
            else:
                self.workflow_sources.append(os.path.abspath(wfs))

        log = "sources="
        for f in self.workflow_sources:
            log += f
        logger.debug(log)
        self.snakefile = workflow.main_snakefile
        self.envvars = {e: os.environ[e] for e in workflow.envvars}
        if self.envvars:
            logger.debug("envvars = %s" % str(self.envvars))
        self.tibanna_sfn = tibanna_sfn
        if precommand:
            self.precommand = precommand
        else:
            self.precommand = ""
        self.s3_bucket = workflow.default_remote_prefix.split("/")[0]
        self.s3_subdir = re.sub(
            "^{}/".format(self.s3_bucket), "", workflow.default_remote_prefix
        )
        logger.debug("precommand= " + self.precommand)
        logger.debug("bucket=" + self.s3_bucket)
        logger.debug("subdir=" + self.s3_subdir)
        self.quiet = quiet

        super().__init__(
            workflow,
            dag,
            cores,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            local_input=local_input,
            restart_times=restart_times,
            assume_shared_fs=False,
            max_status_checks_per_second=max_status_checks_per_second,
            disable_default_remote_provider_args=True,
            disable_default_resources_args=True,
            disable_envvar_declarations=True,
        )
        self.container_image = container_image or get_container_image()
        logger.info(f"Using {self.container_image} for Tibanna jobs.")
        self.tibanna_config = tibanna_config

    def shutdown(self):
        # perform additional steps on shutdown if necessary
        logger.debug("shutting down Tibanna executor")
        super().shutdown()

    def cancel(self):
        from tibanna.core import API

        for j in self.active_jobs:
            logger.info("killing job {}".format(j.jobname))
            while True:
                try:
                    res = API().kill(j.exec_arn)
                    if not self.quiet:
                        print(res)
                    break
                except KeyboardInterrupt:
                    pass
        self.shutdown()

    def split_filename(self, filename, checkdir=None):
        f = os.path.abspath(filename)
        if checkdir:
            checkdir = checkdir.rstrip("/")
            if f.startswith(checkdir):
                fname = re.sub("^{}/".format(checkdir), "", f)
                fdir = checkdir
            else:
                direrrmsg = (
                    "All source files including Snakefile, "
                    + "conda env files, and rule script files "
                    + "must be in the same working directory: {} vs {}"
                )
                raise WorkflowError(direrrmsg.format(checkdir, f))
        else:
            fdir, fname = os.path.split(f)
        return fname, fdir

    def remove_prefix(self, s):
        return re.sub("^{}/{}/".format(self.s3_bucket, self.s3_subdir), "", s)

    def get_job_targets(self, job):
        def handle_target(target):
            if isinstance(target, _IOFile) and target.remote_object.provider.is_default:
                return self.remove_prefix(target)
            else:
                return target

        return [handle_target(target) for target in job.get_targets()]

    def get_snakefile(self):
        return os.path.basename(self.snakefile)

    def add_command(self, job, tibanna_args, tibanna_config):
        # format command
        command = self.format_job_exec(job)

        if self.precommand:
            command = self.precommand + "; " + command
        logger.debug("command = " + str(command))
        tibanna_args.command = command

    def add_workflow_files(self, job, tibanna_args):
        snakefile_fname, snakemake_dir = self.split_filename(self.snakefile)
        snakemake_child_fnames = []
        for src in self.workflow_sources:
            src_fname, _ = self.split_filename(src, snakemake_dir)
            if src_fname != snakefile_fname:  # redundant
                snakemake_child_fnames.append(src_fname)
        # change path for config files
        self.workflow.overwrite_configfiles = [
            self.split_filename(cf, snakemake_dir)[0]
            for cf in self.workflow.overwrite_configfiles
        ]
        tibanna_args.snakemake_directory_local = snakemake_dir
        tibanna_args.snakemake_main_filename = snakefile_fname
        tibanna_args.snakemake_child_filenames = list(set(snakemake_child_fnames))

    def adjust_filepath(self, f):
        if not hasattr(f, "remote_object"):
            rel = self.remove_prefix(f)  # log/benchmark
        elif (
            hasattr(f.remote_object, "provider") and f.remote_object.provider.is_default
        ):
            rel = self.remove_prefix(f)
        else:
            rel = f
        return rel

    def make_tibanna_input(self, job):
        from tibanna import ec2_utils, core as tibanna_core

        # input & output
        # Local snakemake command here must be run with --default-remote-prefix
        # and --default-remote-provider (forced) but on VM these options will be removed.
        # The snakemake on the VM will consider these input and output as not remote.
        # They files are transferred to the container by Tibanna before running snakemake.
        # In short, the paths on VM must be consistent with what's in Snakefile.
        # but the actual location of the files is on the S3 bucket/prefix.
        # This mapping info must be passed to Tibanna.
        for i in job.input:
            logger.debug("job input " + str(i))
            logger.debug("job input is remote= " + ("true" if i.is_remote else "false"))
            if hasattr(i.remote_object, "provider"):
                logger.debug(
                    " is remote default= "
                    + ("true" if i.remote_object.provider.is_default else "false")
                )
        for o in job.expanded_output:
            logger.debug("job output " + str(o))
            logger.debug(
                "job output is remote= " + ("true" if o.is_remote else "false")
            )
            if hasattr(o.remote_object, "provider"):
                logger.debug(
                    " is remote default= "
                    + ("true" if o.remote_object.provider.is_default else "false")
                )
        file_prefix = (
            "file:///data1/snakemake"  # working dir inside snakemake container on VM
        )
        input_source = dict()
        for ip in job.input:
            ip_rel = self.adjust_filepath(ip)
            input_source[os.path.join(file_prefix, ip_rel)] = "s3://" + ip
        output_target = dict()
        output_all = [eo for eo in job.expanded_output]
        if job.log:
            if isinstance(job.log, list):
                output_all.extend([str(_) for _ in job.log])
            else:
                output_all.append(str(job.log))
        if hasattr(job, "benchmark") and job.benchmark:
            if isinstance(job.benchmark, list):
                output_all.extend([str(_) for _ in job.benchmark])
            else:
                output_all.append(str(job.benchmark))
        for op in output_all:
            op_rel = self.adjust_filepath(op)
            output_target[os.path.join(file_prefix, op_rel)] = "s3://" + op

        # mem & cpu
        mem = job.resources["mem_mb"] / 1024 if "mem_mb" in job.resources.keys() else 1
        cpu = job.threads

        # jobid, grouping, run_name
        jobid = tibanna_core.create_jobid()
        if job.is_group():
            run_name = "snakemake-job-%s-group-%s" % (str(jobid), str(job.groupid))
        else:
            run_name = "snakemake-job-%s-rule-%s" % (str(jobid), str(job.rule))

        # tibanna input
        tibanna_config = {
            "run_name": run_name,
            "mem": mem,
            "cpu": cpu,
            "ebs_size": math.ceil(job.resources["disk_mb"] / 1024),
            "log_bucket": self.s3_bucket,
        }
        logger.debug("additional tibanna config: " + str(self.tibanna_config))
        if self.tibanna_config:
            tibanna_config.update(self.tibanna_config)
        tibanna_args = ec2_utils.Args(
            output_S3_bucket=self.s3_bucket,
            language="snakemake",
            container_image=self.container_image,
            input_files=input_source,
            output_target=output_target,
            input_env=self.envvars,
        )
        self.add_workflow_files(job, tibanna_args)
        self.add_command(job, tibanna_args, tibanna_config)
        tibanna_input = {
            "jobid": jobid,
            "config": tibanna_config,
            "args": tibanna_args.as_dict(),
        }
        logger.debug(json.dumps(tibanna_input, indent=4))
        return tibanna_input

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        logger.info("running job using Tibanna...")
        from tibanna.core import API

        super()._run(job)

        # submit job here, and obtain job ids from the backend
        tibanna_input = self.make_tibanna_input(job)
        jobid = tibanna_input["jobid"]
        exec_info = API().run_workflow(
            tibanna_input,
            sfn=self.tibanna_sfn,
            verbose=not self.quiet,
            jobid=jobid,
            open_browser=False,
            sleep=0,
        )
        exec_arn = exec_info.get("_tibanna", {}).get("exec_arn", "")
        jobname = tibanna_input["config"]["run_name"]
        jobid = tibanna_input["jobid"]

        # register job as active, using your own namedtuple.
        # The namedtuple must at least contain the attributes
        # job, jobid, callback, error_callback.
        self.active_jobs.append(
            TibannaJob(job, jobname, jobid, exec_arn, callback, error_callback)
        )

    def _wait_for_jobs(self):
        # busy wait on job completion
        # This is only needed if your backend does not allow to use callbacks
        # for obtaining job status.
        from tibanna.core import API

        while True:
            # always use self.lock to avoid race conditions
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for j in active_jobs:
                # use self.status_rate_limiter to avoid too many API calls.
                with self.status_rate_limiter:
                    if j.exec_arn:
                        status = API().check_status(j.exec_arn)
                    else:
                        status = "FAILED_AT_SUBMISSION"
                    if not self.quiet or status != "RUNNING":
                        logger.debug("job %s: %s" % (j.jobname, status))
                    if status == "RUNNING":
                        still_running.append(j)
                    elif status == "SUCCEEDED":
                        j.callback(j.job)
                    else:
                        j.error_callback(j.job)
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


def run_wrapper(
    job_rule,
    input,
    output,
    params,
    wildcards,
    threads,
    resources,
    log,
    benchmark,
    benchmark_repeats,
    conda_env,
    container_img,
    singularity_args,
    env_modules,
    use_singularity,
    linemaps,
    debug,
    cleanup_scripts,
    shadow_dir,
    jobid,
    edit_notebook,
    conda_base_path,
    basedir,
    runtime_sourcecache_path,
):
    """
    Wrapper around the run method that handles exceptions and benchmarking.

    Arguments
    job_rule   -- the ``job.rule`` member
    input      -- a list of input files
    output     -- a list of output files
    wildcards  -- so far processed wildcards
    threads    -- usable threads
    log        -- a list of log files
    shadow_dir -- optional shadow directory root
    """
    # get shortcuts to job_rule members
    run = job_rule.run_func
    version = job_rule.version
    rule = job_rule.name
    is_shell = job_rule.shellcmd is not None

    if os.name == "posix" and debug:
        sys.stdin = open("/dev/stdin")

    if benchmark is not None:
        from snakemake.benchmark import (
            BenchmarkRecord,
            benchmarked,
            write_benchmark_records,
        )

    # Change workdir if shadow defined and not using singularity.
    # Otherwise, we do the change from inside the container.
    passed_shadow_dir = None
    if use_singularity and container_img:
        passed_shadow_dir = shadow_dir
        shadow_dir = None

    try:
        with change_working_directory(shadow_dir):
            if benchmark:
                bench_records = []
                for bench_iteration in range(benchmark_repeats):
                    # Determine whether to benchmark this process or do not
                    # benchmarking at all.  We benchmark this process unless the
                    # execution is done through the ``shell:``, ``script:``, or
                    # ``wrapper:`` stanza.
                    is_sub = (
                        job_rule.shellcmd
                        or job_rule.script
                        or job_rule.wrapper
                        or job_rule.cwl
                    )
                    if is_sub:
                        # The benchmarking through ``benchmarked()`` is started
                        # in the execution of the shell fragment, script, wrapper
                        # etc, as the child PID is available there.
                        bench_record = BenchmarkRecord()
                        run(
                            input,
                            output,
                            params,
                            wildcards,
                            threads,
                            resources,
                            log,
                            version,
                            rule,
                            conda_env,
                            container_img,
                            singularity_args,
                            use_singularity,
                            env_modules,
                            bench_record,
                            jobid,
                            is_shell,
                            bench_iteration,
                            cleanup_scripts,
                            passed_shadow_dir,
                            edit_notebook,
                            conda_base_path,
                            basedir,
                            runtime_sourcecache_path,
                        )
                    else:
                        # The benchmarking is started here as we have a run section
                        # and the generated Python function is executed in this
                        # process' thread.
                        with benchmarked() as bench_record:
                            run(
                                input,
                                output,
                                params,
                                wildcards,
                                threads,
                                resources,
                                log,
                                version,
                                rule,
                                conda_env,
                                container_img,
                                singularity_args,
                                use_singularity,
                                env_modules,
                                bench_record,
                                jobid,
                                is_shell,
                                bench_iteration,
                                cleanup_scripts,
                                passed_shadow_dir,
                                edit_notebook,
                                conda_base_path,
                                basedir,
                                runtime_sourcecache_path,
                            )
                    # Store benchmark record for this iteration
                    bench_records.append(bench_record)
            else:
                run(
                    input,
                    output,
                    params,
                    wildcards,
                    threads,
                    resources,
                    log,
                    version,
                    rule,
                    conda_env,
                    container_img,
                    singularity_args,
                    use_singularity,
                    env_modules,
                    None,
                    jobid,
                    is_shell,
                    None,
                    cleanup_scripts,
                    passed_shadow_dir,
                    edit_notebook,
                    conda_base_path,
                    basedir,
                    runtime_sourcecache_path,
                )
    except (KeyboardInterrupt, SystemExit) as e:
        # Re-raise the keyboard interrupt in order to record an error in the
        # scheduler but ignore it
        raise e
    except (Exception, BaseException) as ex:
        # this ensures that exception can be re-raised in the parent thread
        origin = get_exception_origin(ex, linemaps)
        if origin is not None:
            log_verbose_traceback(ex)
            lineno, file = origin
            raise RuleException(
                format_error(
                    ex, lineno, linemaps=linemaps, snakefile=file, show_traceback=True
                )
            )
        else:
            # some internal bug, just reraise
            raise ex

    if benchmark is not None:
        try:
            write_benchmark_records(bench_records, benchmark)
        except (Exception, BaseException) as ex:
            raise WorkflowError(ex)
