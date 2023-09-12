__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import ABC, abstractmethod
import asyncio
import os
import sys
import contextlib
import time
import json
import stat
import shutil
import shlex
import threading
import concurrent.futures
import subprocess
import tempfile
from functools import partial
from collections import namedtuple
import base64
from typing import List
import uuid
import re
import math

from snakemake_interface_executor_plugins.executors.base import AbstractExecutor
from snakemake_interface_executor_plugins.executors.real import RealExecutor
from snakemake_interface_executor_plugins.executors.remote import RemoteExecutor
from snakemake_interface_executor_plugins.dag import DAGExecutorInterface
from snakemake_interface_executor_plugins.workflow import WorkflowExecutorInterface
from snakemake_interface_executor_plugins.persistence import StatsExecutorInterface
from snakemake_interface_executor_plugins.logging import LoggerExecutorInterface
from snakemake_interface_executor_plugins.jobs import (
    ExecutorJobInterface,
    SingleJobExecutorInterface,
    GroupJobExecutorInterface,
)
from snakemake_interface_executor_plugins.utils import sleep
from snakemake_interface_executor_plugins.utils import ExecMode

from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.stats import Stats
from snakemake.utils import makedirs
from snakemake.io import get_wildcard_names, Wildcards
from snakemake.exceptions import print_exception, get_exception_origin
from snakemake.exceptions import format_error, RuleException, log_verbose_traceback
from snakemake.exceptions import (
    WorkflowError,
    SpawnedJobError,
    CacheMissException,
)
from snakemake.common import (
    get_container_image,
    get_uuid,
    async_lock,
)


class DryrunExecutor(AbstractExecutor):
    def get_exec_mode(self):
        raise NotImplementedError()

    def printjob(self, job: ExecutorJobInterface):
        super().printjob(job)
        if job.is_group():
            for j in job.jobs:
                self.printcache(j)
        else:
            self.printcache(job)

    def printcache(self, job: ExecutorJobInterface):
        cache_mode = self.workflow.get_cache_mode(job.rule)
        if cache_mode:
            if self.workflow.output_file_cache.exists(job, cache_mode):
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

    def cancel(self):
        pass

    def shutdown(self):
        pass

    def handle_job_success(self, job: ExecutorJobInterface):
        pass

    def handle_job_error(self, job: ExecutorJobInterface):
        pass


class TouchExecutor(RealExecutor):
    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
    ):
        super().__init__(
            workflow,
            dag,
            stats,
            logger,
            executor_settings=None,
        )

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
        super()._run(job)
        try:
            # Touching of output files will be done by handle_job_success
            time.sleep(0.1)
            callback(job)
        except OSError as ex:
            print_exception(ex, self.workflow.linemaps)
            error_callback(job)

    def get_exec_mode(self):
        raise NotImplementedError()

    def handle_job_success(self, job: ExecutorJobInterface):
        super().handle_job_success(job, ignore_missing_output=True)

    def cancel(self):
        pass

    def shutdown(self):
        pass

    def get_python_executable(self):
        raise NotImplementedError()


_ProcessPoolExceptions = (KeyboardInterrupt,)
try:
    from concurrent.futures.process import BrokenProcessPool

    _ProcessPoolExceptions = (KeyboardInterrupt, BrokenProcessPool)
except ImportError:
    pass


class CPUExecutor(RealExecutor):
    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
        cores: int,
        use_threads=False,
    ):
        super().__init__(
            workflow,
            dag,
            stats,
            logger,
            executor_settings=None,
            job_core_limit=cores,
        )

        self.use_threads = use_threads

        # Zero thread jobs do not need a thread, but they occupy additional workers.
        # Hence we need to reserve additional workers for them.
        workers = cores + 5 if cores is not None else 5
        self.workers = workers
        self.pool = concurrent.futures.ThreadPoolExecutor(max_workers=self.workers)

    def get_exec_mode(self):
        return ExecMode.subprocess

    @property
    def job_specific_local_groupid(self):
        return False

    def get_job_exec_prefix(self, job: ExecutorJobInterface):
        return f"cd {shlex.quote(self.workflow.workdir_init)}"

    def get_python_executable(self):
        return sys.executable

    def get_envvar_declarations(self):
        return ""

    def get_job_args(self, job: ExecutorJobInterface, **kwargs):
        return f"{super().get_job_args(job, **kwargs)} --quiet"

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
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

    def job_args_and_prepare(self, job: ExecutorJobInterface):
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

    def run_single_job(self, job: SingleJobExecutorInterface):
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

    def run_group_job(self, job: GroupJobExecutorInterface):
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

    def spawn_job(self, job: SingleJobExecutorInterface):
        cmd = self.format_job_exec(job)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            raise SpawnedJobError()

    def cached_or_run(self, job: SingleJobExecutorInterface, run_func, *args):
        """
        Either retrieve result from cache, or run job with given function.
        """
        cache_mode = self.workflow.get_cache_mode(job.rule)
        try:
            if cache_mode:
                self.workflow.output_file_cache.fetch(job, cache_mode)
                return
        except CacheMissException:
            pass
        run_func(*args)
        if cache_mode:
            self.workflow.output_file_cache.store(job, cache_mode)

    def shutdown(self):
        self.pool.shutdown()

    def cancel(self):
        self.pool.shutdown()

    def _callback(
        self, job: SingleJobExecutorInterface, callback, error_callback, future
    ):
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
        except BaseException as ex:
            self.print_job_error(job)
            if self.workflow.verbose or (not job.is_group() and not job.is_shell):
                print_exception(ex, self.workflow.linemaps)
            error_callback(job)

    def handle_job_success(self, job: ExecutorJobInterface):
        super().handle_job_success(job)

    def handle_job_error(self, job: ExecutorJobInterface):
        super().handle_job_error(job)
        if not self.keepincomplete:
            job.cleanup()
            self.workflow.persistence.cleanup(job)


GenericClusterJob = namedtuple(
    "GenericClusterJob",
    "job jobid callback error_callback jobscript jobfinished jobfailed",
)


class GenericClusterExecutor(RemoteExecutor):
    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
        submitcmd="qsub",
        statuscmd=None,
        cancelcmd=None,
        cancelnargs=None,
        sidecarcmd=None,
        jobname="snakejob.{rulename}.{jobid}.sh",
        max_status_checks_per_second=1,
    ):
        self.submitcmd = submitcmd
        if not workflow.assume_shared_fs and statuscmd is None:
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
            stats,
            logger,
            None,
            jobname=jobname,
            max_status_checks_per_second=max_status_checks_per_second,
        )

        self.sidecar_vars = None
        if self.sidecarcmd:
            self._launch_sidecar()

        if not statuscmd and not self.assume_shared_fs:
            raise WorkflowError(
                "If no shared filesystem is used, you have to "
                "specify a cluster status command."
            )

    def get_job_exec_prefix(self, job: ExecutorJobInterface):
        if self.assume_shared_fs:
            return f"cd {shlex.quote(self.workflow.workdir_init)}"
        else:
            return ""

    def get_job_exec_suffix(self, job: ExecutorJobInterface):
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

    def get_jobfinished_marker(self, job: ExecutorJobInterface):
        return os.path.join(self.tmpdir, f"{job.jobid}.jobfinished")

    def get_jobfailed_marker(self, job: ExecutorJobInterface):
        return os.path.join(self.tmpdir, f"{job.jobid}.jobfailed")

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

    def register_job(self, job: ExecutorJobInterface):
        # Do not register job here.
        # Instead do it manually once the jobid is known.
        pass

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
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
            submitcmd = job.format_wildcards(self.submitcmd, dependencies=deps)
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

    async def _wait_for_jobs(self):
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
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            # logger.debug("Checking status of {} jobs.".format(len(active_jobs)))
            for active_job in active_jobs:
                async with self.status_rate_limiter:
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
            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            await sleep()


SynchronousClusterJob = namedtuple(
    "SynchronousClusterJob", "job jobid callback error_callback jobscript process"
)


class SynchronousClusterExecutor(RemoteExecutor):
    """
    invocations like "qsub -sync y" (SGE) or "bsub -K" (LSF) are
    synchronous, blocking the foreground thread and returning the
    remote exit code at remote exit.
    """

    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
        submitcmd="qsub",
        jobname="snakejob.{rulename}.{jobid}.sh",
    ):
        super().__init__(
            workflow,
            dag,
            stats,
            logger,
            None,
            jobname=jobname,
            max_status_checks_per_second=10,
        )
        self.submitcmd = submitcmd
        self.external_jobid = dict()

    def get_job_exec_prefix(self, job):
        if self.assume_shared_fs:
            return f"cd {shlex.quote(self.workflow.workdir_init)}"
        else:
            return ""

    def cancel(self):
        logger.info("Will exit after finishing currently running jobs.")
        self.shutdown()

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
        super()._run(job)

        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        deps = " ".join(
            self.external_jobid[f] for f in job.input if f in self.external_jobid
        )
        try:
            submitcmd = job.format_wildcards(self.submitcmd, dependencies=deps)
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

    async def _wait_for_jobs(self):
        while True:
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for active_job in active_jobs:
                async with self.status_rate_limiter:
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
            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            await sleep()


DRMAAClusterJob = namedtuple(
    "DRMAAClusterJob", "job jobid callback error_callback jobscript"
)


class DRMAAExecutor(RemoteExecutor):
    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
        jobname="snakejob.{rulename}.{jobid}.sh",
        drmaa_args="",
        drmaa_log_dir=None,
        max_status_checks_per_second=1,
    ):
        super().__init__(
            workflow,
            dag,
            stats,
            logger,
            None,
            jobname=jobname,
            max_status_checks_per_second=max_status_checks_per_second,
        )
        try:
            import drmaa
        except ImportError:
            raise WorkflowError(
                "Python support for DRMAA is not installed. "
                "Please install it, e.g. with easy_install3 --user drmaa"
            )
        except RuntimeError as e:
            raise WorkflowError(f"Error loading drmaa support:\n{e}")
        self.session = drmaa.Session()
        self.drmaa_args = drmaa_args
        self.drmaa_log_dir = drmaa_log_dir
        self.session.initialize()
        self.submitted = list()

    def get_job_exec_prefix(self, job: ExecutorJobInterface):
        if self.assume_shared_fs:
            return f"cd {shlex.quote(self.workflow.workdir_init)}"
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

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
        super()._run(job)
        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        try:
            drmaa_args = job.format_wildcards(self.drmaa_args)
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
            print_exception(WorkflowError(f"DRMAA Error: {e}"), self.workflow.linemaps)
            error_callback(job)
            return
        logger.info(f"Submitted DRMAA job {job.jobid} with external jobid {jobid}.")
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

    async def _wait_for_jobs(self):
        import drmaa

        suspended_msg = set()

        while True:
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for active_job in active_jobs:
                async with self.status_rate_limiter:
                    try:
                        retval = self.session.jobStatus(active_job.jobid)
                    except drmaa.ExitTimeoutException as e:
                        # job still active
                        still_running.append(active_job)
                        continue
                    except (drmaa.InternalException, Exception) as e:
                        print_exception(
                            WorkflowError(f"DRMAA Error: {e}"),
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

            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            await sleep()


@contextlib.contextmanager
def change_working_directory(directory=None):
    """Change working directory in execution context if provided."""
    if directory:
        try:
            saved_directory = os.getcwd()
            logger.info(f"Changing to shadow directory: {directory}")
            os.chdir(directory)
            yield
        finally:
            os.chdir(saved_directory)
    else:
        yield


KubernetesJob = namedtuple(
    "KubernetesJob", "job jobid callback error_callback kubejob jobscript"
)


class KubernetesExecutor(RemoteExecutor):
    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
        namespace,
        container_image=None,
        k8s_cpu_scalar=1.0,
        k8s_service_account_name=None,
        jobname="{rulename}.{jobid}",
    ):
        self.workflow = workflow

        super().__init__(
            workflow,
            dag,
            stats,
            logger,
            None,
            jobname=jobname,
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

        self.k8s_cpu_scalar = k8s_cpu_scalar
        self.k8s_service_account_name = k8s_service_account_name
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

    def get_job_exec_prefix(self, job: ExecutorJobInterface):
        return "cp -rf /source/. ."

    def register_secret(self):
        import kubernetes.client

        secret = kubernetes.client.V1Secret()
        secret.metadata = kubernetes.client.V1ObjectMeta()
        # create a random uuid
        secret.metadata.name = self.run_namespace
        secret.type = "Opaque"
        secret.data = {}
        for i, f in enumerate(self.dag.get_sources()):
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
                key = f"f{i}"

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
                logger.warning(f"  * File: {k}, original size: {v}")

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

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
        import kubernetes.client

        super()._run(job)
        exec_job = self.format_job_exec(job)

        # Kubernetes silently does not submit a job if the name is too long
        # therefore, we ensure that it is not longer than snakejob+uuid.
        jobid = "snakejob-{}".format(
            get_uuid(f"{self.run_namespace}-{job.jobid}-{job.attempt}")
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
            kubernetes.client.V1VolumeMount(name="workdir", mount_path="/workdir"),
            kubernetes.client.V1VolumeMount(name="source", mount_path="/source"),
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
        # Add service account name if provided
        if self.k8s_service_account_name:
            body.spec.service_account_name = self.k8s_service_account_name

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
        logger.debug(f"job resources:  {dict(job.resources)}")
        container.resources = kubernetes.client.V1ResourceRequirements()
        container.resources.requests = {}
        container.resources.requests["cpu"] = "{}m".format(
            int(job.resources["_cores"] * self.k8s_cpu_scalar * 1000)
        )
        if "mem_mb" in job.resources.keys():
            container.resources.requests["memory"] = "{}M".format(
                job.resources["mem_mb"]
            )
        if "disk_mb" in job.resources.keys():
            disk_mb = int(job.resources.get("disk_mb", 1024))
            container.resources.requests["ephemeral-storage"] = f"{disk_mb}M"

        logger.debug(f"k8s pod resources: {container.resources.requests}")

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

    async def _wait_for_jobs(self):
        import kubernetes

        while True:
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for j in active_jobs:
                async with self.status_rate_limiter:
                    logger.debug(f"Checking status for pod {j.jobid}")
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
            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            await sleep()


TibannaJob = namedtuple(
    "TibannaJob", "job jobname jobid exec_arn callback error_callback"
)


class TibannaExecutor(RemoteExecutor):
    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
        tibanna_sfn,
        precommand="",
        tibanna_config=False,
        container_image=None,
        max_status_checks_per_second=1,
    ):
        super().__init__(
            workflow,
            dag,
            stats,
            logger,
            None,
            max_status_checks_per_second=max_status_checks_per_second,
            disable_default_remote_provider_args=True,
            disable_default_resources_args=True,
            disable_envvar_declarations=True,
        )
        self.workflow = workflow
        self.workflow_sources = []
        for wfs in dag.get_sources():
            if os.path.isdir(wfs):
                for dirpath, dirnames, filenames in os.walk(wfs):
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
            f"^{self.s3_bucket}/", "", workflow.default_remote_prefix
        )
        logger.debug("precommand= " + self.precommand)
        logger.debug("bucket=" + self.s3_bucket)
        logger.debug("subdir=" + self.s3_subdir)
        self.quiet = workflow.quiet

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
            logger.info(f"killing job {j.jobname}")
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
                fname = re.sub(f"^{checkdir}/", "", f)
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
        return re.sub(f"^{self.s3_bucket}/{self.s3_subdir}/", "", s)

    def get_snakefile(self):
        return os.path.basename(self.snakefile)

    def add_command(self, job: ExecutorJobInterface, tibanna_args, tibanna_config):
        # format command
        command = self.format_job_exec(job)

        if self.precommand:
            command = self.precommand + "; " + command
        logger.debug("command = " + str(command))
        tibanna_args.command = command

    def add_workflow_files(self, job: ExecutorJobInterface, tibanna_args):
        snakefile_fname, snakemake_dir = self.split_filename(self.snakefile)
        snakemake_child_fnames = []
        for src in self.workflow_sources:
            src_fname, _ = self.split_filename(src, snakemake_dir)
            if src_fname != snakefile_fname:  # redundant
                snakemake_child_fnames.append(src_fname)
        # change path for config files
        # TODO - this is a hacky way to do this
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

    def make_tibanna_input(self, job: ExecutorJobInterface):
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
            run_name = f"snakemake-job-{str(jobid)}-group-{str(job.groupid)}"
        else:
            run_name = f"snakemake-job-{str(jobid)}-rule-{str(job.rule)}"

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

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
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

    async def _wait_for_jobs(self):
        # busy wait on job completion
        # This is only needed if your backend does not allow to use callbacks
        # for obtaining job status.
        from tibanna.core import API

        while True:
            # always use self.lock to avoid race conditions
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for j in active_jobs:
                # use self.status_rate_limiter to avoid too many API calls.
                async with self.status_rate_limiter:
                    if j.exec_arn:
                        status = API().check_status(j.exec_arn)
                    else:
                        status = "FAILED_AT_SUBMISSION"
                    if not self.quiet or status != "RUNNING":
                        logger.debug(f"job {j.jobname}: {status}")
                    if status == "RUNNING":
                        still_running.append(j)
                    elif status == "SUCCEEDED":
                        j.callback(j.job)
                    else:
                        j.error_callback(j.job)
            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            await sleep()


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
    except BaseException as ex:
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
        except BaseException as ex:
            raise WorkflowError(ex)
