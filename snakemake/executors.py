__author__ = "Johannes Köster"
__contributors__ = ["David Alexander"]
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

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
from functools import partial
from itertools import chain
from collections import namedtuple
from tempfile import mkdtemp
import random
import base64
import uuid

from snakemake.jobs import Job
from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.stats import Stats
from snakemake.utils import format, Unformattable, makedirs
from snakemake.io import get_wildcard_names, Wildcards
from snakemake.exceptions import print_exception, get_exception_origin
from snakemake.exceptions import format_error, RuleException, log_verbose_traceback
from snakemake.exceptions import ClusterJobException, ProtectedOutputException, WorkflowError, ImproperShadowException, SpawnedJobError
from snakemake.common import Mode, __version__, get_container_image, get_uuid


def sleep():
    # do not sleep on CI. In that case we just want to quickly test everything.
    if os.environ.get("CIRCLECI") != "true":
        time.sleep(10)


class AbstractExecutor:
    def __init__(self, workflow, dag,
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 printthreads=True,
                 latency_wait=3):
        self.workflow = workflow
        self.dag = dag
        self.quiet = quiet
        self.printreason = printreason
        self.printshellcmds = printshellcmds
        self.printthreads = printthreads
        self.latency_wait = latency_wait

    def get_default_remote_provider_args(self):
        if self.workflow.default_remote_provider:
            return (
                " --default-remote-provider {} "
                "--default-remote-prefix {} ").format(
                    self.workflow.default_remote_provider.__module__.split(".")[-1],
                    self.workflow.default_remote_prefix)
        return ""

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
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
    pass


class RealExecutor(AbstractExecutor):
    def __init__(self, workflow, dag,
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 assume_shared_fs=True):
        super().__init__(workflow, dag,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait)
        self.assume_shared_fs = assume_shared_fs
        self.stats = Stats()
        self.snakefile = workflow.snakefile

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
                "directory {}".format(e, self.workflow.persistence.path))

    def handle_job_success(self, job,
                           upload_remote=True,
                           handle_log=True,
                           handle_touch=True,
                           ignore_missing_output=False):
        job.postprocess(upload_remote=upload_remote,
                        handle_log=handle_log,
                        handle_touch=handle_touch,
                        ignore_missing_output=ignore_missing_output,
                        latency_wait=self.latency_wait,
                        assume_shared_fs=self.assume_shared_fs)
        self.stats.report_job_end(job)

    def handle_job_error(self, job, upload_remote=True):
        job.postprocess(error=True,
                        assume_shared_fs=self.assume_shared_fs,
                        latency_wait=self.latency_wait)

    def format_job_pattern(self, pattern, job=None, **kwargs):
        overwrite_workdir = []
        if self.workflow.overwrite_workdir:
            overwrite_workdir.extend(("--directory", self.workflow.overwrite_workdir))

        overwrite_config = []
        if self.workflow.overwrite_configfile:
            overwrite_config.extend(("--configfile", self.workflow.overwrite_configfile))
        if self.workflow.config_args:
            overwrite_config.append("--config")
            overwrite_config.extend(self.workflow.config_args)

        printshellcmds = ""
        if self.workflow.printshellcmds:
            printshellcmds = "-p"

        return format(pattern,
                      job=job,
                      attempt=job.attempt,
                      overwrite_workdir=overwrite_workdir,
                      overwrite_config=overwrite_config,
                      printshellcmds=printshellcmds,
                      workflow=self.workflow,
                      snakefile=self.snakefile,
                      cores=self.cores,
                      benchmark_repeats=job.benchmark_repeats if not job.is_group() else None,
                      target=job.get_targets(),
                      **kwargs)


class TouchExecutor(RealExecutor):
    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        super()._run(job)
        try:
            #Touching of output files will be done by handle_job_success
            time.sleep(0.1)
            callback(job)
        except OSError as ex:
            print_exception(ex, self.workflow.linemaps)
            error_callback(job)

    def handle_job_success(self, job):
        super().handle_job_success(job, ignore_missing_output=True)


_ProcessPoolExceptions = (KeyboardInterrupt, )
try:
    from concurrent.futures.process import BrokenProcessPool
    _ProcessPoolExceptions = (KeyboardInterrupt, BrokenProcessPool)
except ImportError:
    pass


class CPUExecutor(RealExecutor):
    def __init__(self, workflow, dag, workers,
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 use_threads=False,
                 latency_wait=3,
                 cores=1):
        super().__init__(workflow, dag,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait)

        self.exec_job = '\\\n'.join((
            'cd {workflow.workdir_init} && ',
            '{sys.executable} -m snakemake {target} --snakefile {snakefile} ',
            '--force -j{cores} --keep-target-files --keep-remote ',
            '--attempt {attempt} ',
            '--force-use-threads --wrapper-prefix {workflow.wrapper_prefix} ',
            '--latency-wait {latency_wait} ',
            self.get_default_remote_provider_args(),
            '{overwrite_workdir} {overwrite_config} {printshellcmds} ',
            '--notemp --quiet --no-hooks --nolock --mode {} '.format(Mode.subprocess)))

        if self.workflow.shadow_prefix:
            self.exec_job += " --shadow-prefix {} ".format(
                self.workflow.shadow_prefix)
        if self.workflow.use_conda:
            self.exec_job += " --use-conda "
            if self.workflow.conda_prefix:
                self.exec_job += " --conda-prefix {} ".format(
                    self.workflow.conda_prefix)
        if self.workflow.use_singularity:
            self.exec_job += " --use-singularity "
            if self.workflow.singularity_prefix:
                self.exec_job += " --singularity-prefix {} ".format(
                    self.workflow.singularity_prefix)
            if self.workflow.singularity_args:
                self.exec_job += " --singularity-args \"{}\"".format(
                    self.workflow.singularity_args)

        self.use_threads = use_threads
        self.cores = cores
        self.pool = concurrent.futures.ThreadPoolExecutor(max_workers=workers + 1)

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        super()._run(job)

        if job.is_group():
            # the future waits for the entire group job
            future = self.pool.submit(self.run_group_job, job)
        else:
            future = self.run_single_job(job)

        future.add_done_callback(partial(self._callback, job, callback,
                                         error_callback))

    def job_args_and_prepare(self, job):
        job.prepare()

        conda_env = job.conda_env_path
        singularity_img = job.singularity_img_path

        benchmark = None
        benchmark_repeats = job.benchmark_repeats or 1
        if job.benchmark is not None:
            benchmark = str(job.benchmark)
        return (job.rule, job.input.plainstrings(),
                job.output.plainstrings(), job.params, job.wildcards,
                job.threads, job.resources, job.log.plainstrings(), benchmark,
                benchmark_repeats, conda_env, singularity_img,
                self.workflow.singularity_args, self.workflow.use_singularity,
                self.workflow.linemaps, self.workflow.debug,
                job.shadow_dir, job.jobid)

    def run_single_job(self, job):
        if self.use_threads or (not job.is_shadow and not job.is_run):
            future = self.pool.submit(
                run_wrapper, *self.job_args_and_prepare(job))
        else:
            # run directive jobs are spawned into subprocesses
            future = self.pool.submit(self.spawn_job, job)
        return future

    def run_group_job(self, job):
        """Run a pipe group job.

        This lets all items run simultaneously."""
        # we only have to consider pipe groups because in local running mode,
        # these are the only groups that will occur

        futures = [self.run_single_job(j) for j in job]

        while True:
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
                        return
            time.sleep(1)

    def spawn_job(self, job):
        exec_job = self.exec_job
        if not job.is_branched:
            exec_job += " --allowed-rules {}".format(" ".join(job.rules))
        cmd = self.format_job_pattern(exec_job, job=job,
                                      _quote_all=True,
                                      latency_wait=self.latency_wait)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            raise SpawnedJobError()

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
            print_exception(ex, self.workflow.linemaps)
            error_callback(job)

    def handle_job_success(self, job):
        super().handle_job_success(job)

    def handle_job_error(self, job):
        super().handle_job_error(job)
        job.cleanup()
        self.workflow.persistence.cleanup(job)


class ClusterExecutor(RealExecutor):

    default_jobscript = "jobscript.sh"

    def __init__(self, workflow, dag, cores,
                 jobname="snakejob.{name}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 cluster_config=None,
                 local_input=None,
                 restart_times=None,
                 exec_job=None,
                 assume_shared_fs=True,
                 max_status_checks_per_second=1):
        from ratelimiter import RateLimiter

        local_input = local_input or []
        super().__init__(workflow, dag,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         assume_shared_fs=assume_shared_fs)

        if not self.assume_shared_fs:
            # use relative path to Snakefile
            self.snakefile = os.path.relpath(workflow.snakefile)

        jobscript = workflow.jobscript
        if jobscript is None:
            jobscript = os.path.join(os.path.dirname(__file__),
                                     self.default_jobscript)
        try:
            with open(jobscript) as f:
                self.jobscript = f.read()
        except IOError as e:
            raise WorkflowError(e)

        if not "jobid" in get_wildcard_names(jobname):
            raise WorkflowError(
                "Defined jobname (\"{}\") has to contain the wildcard {jobid}.")

        if exec_job is None:
            self.exec_job = '\\\n'.join((
                'cd {workflow.workdir_init} && ' if assume_shared_fs else '',
                '{sys.executable} ' if assume_shared_fs else 'python ',
                '-m snakemake {target} --snakefile {snakefile} ',
                '--force -j{cores} --keep-target-files --keep-remote ',
                '--wait-for-files {wait_for_files} --latency-wait {latency_wait} ',
                ' --attempt {attempt} {use_threads} ',
                '--wrapper-prefix {workflow.wrapper_prefix} ',
                '{overwrite_workdir} {overwrite_config} {printshellcmds} --nocolor ',
                '--notemp --no-hooks --nolock --mode {} '.format(Mode.cluster)))
        else:
            self.exec_job = exec_job

        if self.workflow.shadow_prefix:
            self.exec_job += " --shadow-prefix {} ".format(
                self.workflow.shadow_prefix)
        if self.workflow.use_conda:
            self.exec_job += " --use-conda "
            if self.workflow.conda_prefix:
                self.exec_job += " --conda-prefix {} ".format(
                    self.workflow.conda_prefix)
        if self.workflow.use_singularity:
            self.exec_job += " --use-singularity "
            if self.workflow.singularity_prefix:
                self.exec_job += " --singularity-prefix {} ".format(
                    self.workflow.singularity_prefix)
            if self.workflow.singularity_args:
                self.exec_job += " --singularity-args \"{}\"".format(
                    self.workflow.singularity_args)

        self.exec_job += self.get_default_remote_provider_args()

        if not any(dag.dynamic_output_jobs):
            # disable restiction to target rule in case of dynamic rules!
            self.exec_job += " --allowed-rules {rules} "
        self.jobname = jobname
        self._tmpdir = None
        self.cores = cores if cores else ""
        self.cluster_config = cluster_config if cluster_config else dict()

        self.restart_times = restart_times

        self.active_jobs = list()
        self.lock = threading.Lock()
        self.wait = True
        self.wait_thread = threading.Thread(target=self._wait_for_jobs)
        self.wait_thread.daemon = True
        self.wait_thread.start()

        self.max_status_checks_per_second = max_status_checks_per_second

        self.status_rate_limiter = RateLimiter(
            max_calls=self.max_status_checks_per_second,
            period=1)

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
            self._tmpdir = mkdtemp(dir=".snakemake", prefix="tmp.")
        return os.path.abspath(self._tmpdir)

    def get_jobscript(self, job):
        f = job.format_wildcards(self.jobname,
                             rulename=job.name,
                             name=job.name,
                             jobid=job.jobid,
                             cluster=self.cluster_wildcards(job))

        if os.path.sep in f:
            raise WorkflowError("Path separator ({}) found in job name {}. "
                                "This is not supported.".format(
                                os.path.sep, f))

        return os.path.join(self.tmpdir, f)

    def format_job(self, pattern, job, **kwargs):
        wait_for_files = []
        if self.assume_shared_fs:
            wait_for_files.append(self.tmpdir)
            wait_for_files.extend(job.get_wait_for_files())

        format_p = partial(self.format_job_pattern,
                           job=job,
                           properties=job.properties(
                               cluster=self.cluster_params(job)),
                           latency_wait=self.latency_wait,
                           wait_for_files=wait_for_files,
                           **kwargs)
        try:
            return format_p(pattern)
        except KeyError as e:
            raise WorkflowError(
                "Error formatting jobscript: {} not found\n"
                "Make sure that your custom jobscript is up to date.".format(e))

    def write_jobscript(self, job, jobscript, **kwargs):
        # only force threads if this is not a group job
        # otherwise we want proper process handling
        use_threads = "--force-use-threads" if not job.is_group() else ""
        exec_job = self.format_job(self.exec_job,
                                   job,
                                   _quote_all=True,
                                   rules=job.rules,
                                   use_threads=use_threads,
                                   **kwargs)
        content = self.format_job(self.jobscript,
                                  job,
                                  exec_job=exec_job,
                                  **kwargs)
        logger.debug("Jobscript:\n{}".format(content))
        with open(jobscript, "w") as f:
            print(content, file=f)
        os.chmod(jobscript, os.stat(jobscript).st_mode | stat.S_IXUSR)

    def cluster_params(self, job):
        """Return wildcards object for job from cluster_config."""
        cluster = self.cluster_config.get("__default__", dict()).copy()
        cluster.update(self.cluster_config.get(job.name, dict()))
        # Format values with available parameters from the job.
        for key, value in list(cluster.items()):
            if isinstance(value, str):
                cluster[key] = job.format_wildcards(value)

        return cluster

    def cluster_wildcards(self, job):
        return Wildcards(fromdict=self.cluster_params(job))

    def handle_job_success(self, job):
        super().handle_job_success(job, upload_remote=False,
                                   handle_log=False, handle_touch=False)

    def handle_job_error(self, job):
        # TODO what about removing empty remote dirs?? This cannot be decided
        # on the cluster node.
        super().handle_job_error(job, upload_remote=False)
        logger.debug("Cleanup job metadata.")
        # We have to remove metadata here as well.
        # It will be removed by the CPUExecutor in case of a shared FS,
        # but we might not see the removal due to filesystem latency.
        # By removing it again, we make sure that it is gone on the host FS.
        self.workflow.persistence.cleanup(job)


GenericClusterJob = namedtuple("GenericClusterJob", "job jobid callback error_callback jobscript jobfinished jobfailed")


class GenericClusterExecutor(ClusterExecutor):
    def __init__(self, workflow, dag, cores,
                 submitcmd="qsub",
                 statuscmd=None,
                 cluster_config=None,
                 jobname="snakejob.{rulename}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 restart_times=0,
                 assume_shared_fs=True,
                 max_status_checks_per_second=1):

        self.submitcmd = submitcmd
        if not assume_shared_fs and statuscmd is None:
            raise WorkflowError("When no shared filesystem can be assumed, a "
                "status command must be given.")

        self.statuscmd = statuscmd
        self.external_jobid = dict()

        super().__init__(workflow, dag, cores,
                         jobname=jobname,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         cluster_config=cluster_config,
                         restart_times=restart_times,
                         assume_shared_fs=assume_shared_fs,
                         max_status_checks_per_second=max_status_checks_per_second)

        if statuscmd:
            self.exec_job += ' && exit 0 || exit 1'
        elif assume_shared_fs:
            # TODO wrap with watch and touch {jobrunning}
            # check modification date of {jobrunning} in the wait_for_job method
            self.exec_job += ' && touch "{jobfinished}" || (touch "{jobfailed}"; exit 1)'
        else:
            raise WorkflowError("If no shared filesystem is used, you have to "
                                "specify a cluster status command.")


    def cancel(self):
        logger.info("Will exit after finishing currently running jobs.")
        self.shutdown()

    def register_job(self, job):
        # Do not register job here.
        # Instead do it manually once the jobid is known.
        pass

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        super()._run(job)
        workdir = os.getcwd()
        jobid = job.jobid

        jobscript = self.get_jobscript(job)
        jobfinished = os.path.join(self.tmpdir, "{}.jobfinished".format(jobid))
        jobfailed = os.path.join(self.tmpdir, "{}.jobfailed".format(jobid))
        self.write_jobscript(job, jobscript,
                             jobfinished=jobfinished,
                             jobfailed=jobfailed)

        if self.statuscmd:
            ext_jobid = self.dag.incomplete_external_jobid(job)
            if ext_jobid:
                # Job is incomplete and still running.
                # We simply register it and wait for completion or failure.
                logger.info(
                    "Resuming incomplete job {} with external jobid '{}'.".format(
                    jobid, ext_jobid))
                submit_callback(job)
                with self.lock:
                    self.active_jobs.append(
                        GenericClusterJob(job,
                                          ext_jobid,
                                          callback,
                                          error_callback,
                                          jobscript,
                                          jobfinished,
                                          jobfailed))
                return

        deps = " ".join(self.external_jobid[f] for f in job.input
                        if f in self.external_jobid)
        try:
            submitcmd = job.format_wildcards(
                self.submitcmd,
                dependencies=deps,
                cluster=self.cluster_wildcards(job))
        except AttributeError as e:
            raise WorkflowError(str(e),
                                rule=job.rule if not job.is_group() else None)

        try:
            ext_jobid = subprocess.check_output(
                '{submitcmd} "{jobscript}"'.format(submitcmd=submitcmd,
                                                   jobscript=jobscript),
                shell=True).decode().split("\n")
        except subprocess.CalledProcessError as ex:
            logger.error("Error submitting jobscript (exit code {}):\n{}".format(
                    ex.returncode, ex.output.decode()))
            error_callback(job)
            return
        if ext_jobid and ext_jobid[0]:
            ext_jobid = ext_jobid[0]
            self.external_jobid.update((f, ext_jobid) for f in job.output)
            logger.info("Submitted {} {} with external jobid '{}'.".format(
                "group job" if job.is_group() else "job",
                jobid, ext_jobid))
            self.workflow.persistence.started(
                job, external_jobid=ext_jobid)

        submit_callback(job)

        with self.lock:
            self.active_jobs.append(GenericClusterJob(
                job, ext_jobid, callback, error_callback, jobscript,
                jobfinished, jobfailed))

    def _wait_for_jobs(self):
        success = "success"
        failed = "failed"
        running = "running"
        if self.statuscmd is not None:
            def job_status(job):
                try:
                    # this command shall return "success", "failed" or "running"
                    return subprocess.check_output(
                        '{statuscmd} {jobid}'.format(jobid=job.jobid,
                                                     statuscmd=self.statuscmd),
                        shell=True).decode().split("\n")[0]
                except subprocess.CalledProcessError as e:
                    if e.returncode < 0:
                        # Ignore SIGINT and all other issues due to signals
                        # because it will be caused by hitting e.g.
                        # Ctrl-C on the main process or sending killall to
                        # snakemake.
                        # Snakemake will handle the signal in
                        # the master process.
                        pass
                    else:
                        raise WorkflowError("Failed to obtain job status. "
                                            "See above for error message.")
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
            logger.debug("Checking status of {} jobs.".format(len(active_jobs)))
            for active_job in active_jobs:
                with self.status_rate_limiter:
                    status = job_status(active_job)

                    if status == success:
                        active_job.callback(active_job.job)
                    elif status == failed:
                        self.print_job_error(
                            active_job.job,
                            cluster_jobid=active_job.jobid if active_job.jobid else "unknown",
                        )
                        active_job.error_callback(active_job.job)
                    else:
                        still_running.append(active_job)
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


SynchronousClusterJob = namedtuple("SynchronousClusterJob", "job jobid callback error_callback jobscript process")


class SynchronousClusterExecutor(ClusterExecutor):
    """
    invocations like "qsub -sync y" (SGE) or "bsub -K" (LSF) are
    synchronous, blocking the foreground thread and returning the
    remote exit code at remote exit.
    """

    def __init__(self, workflow, dag, cores,
                 submitcmd="qsub",
                 cluster_config=None,
                 jobname="snakejob.{rulename}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 restart_times=0,
                 assume_shared_fs=True):
        super().__init__(workflow, dag, cores,
                         jobname=jobname,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         cluster_config=cluster_config,
                         restart_times=restart_times,
                         assume_shared_fs=assume_shared_fs,
                         max_status_checks_per_second=10)
        self.submitcmd = submitcmd
        self.external_jobid = dict()

    def cancel(self):
        logger.info("Will exit after finishing currently running jobs.")
        self.shutdown()

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        super()._run(job)
        workdir = os.getcwd()
        jobid = job.jobid

        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        deps = " ".join(self.external_jobid[f] for f in job.input
                        if f in self.external_jobid)
        try:
            submitcmd = job.format_wildcards(
                self.submitcmd,
                dependencies=deps,
                cluster=self.cluster_wildcards(job))
        except AttributeError as e:
            raise WorkflowError(str(e),
                                rule=job.rule if not job.is_group() else None)

        process = subprocess.Popen('{submitcmd} "{jobscript}"'.format(submitcmd=submitcmd,
                                           jobscript=jobscript), shell=True)
        submit_callback(job)

        with self.lock:
            self.active_jobs.append(SynchronousClusterJob(
                job, process.pid, callback, error_callback, jobscript,
                process))

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
                        print_exception(
                            ClusterJobException(
                                active_job, self.dag.jobid(active_job.job)),
                            self.workflow.linemaps)
                        active_job.error_callback(active_job.job)
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


DRMAAClusterJob = namedtuple("DRMAAClusterJob", "job jobid callback error_callback jobscript")


class DRMAAExecutor(ClusterExecutor):
    def __init__(self, workflow, dag, cores,
                 jobname="snakejob.{rulename}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 drmaa_args="",
                 drmaa_log_dir=None,
                 latency_wait=3,
                 cluster_config=None,
                 restart_times=0,
                 assume_shared_fs=True,
                 max_status_checks_per_second=1):
        super().__init__(workflow, dag, cores,
                         jobname=jobname,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         cluster_config=cluster_config,
                         restart_times=restart_times,
                         assume_shared_fs=assume_shared_fs,
                         max_status_checks_per_second=max_status_checks_per_second)
        try:
            import drmaa
        except ImportError:
            raise WorkflowError(
                "Python support for DRMAA is not installed. "
                "Please install it, e.g. with easy_install3 --user drmaa")
        except RuntimeError as e:
            raise WorkflowError("Error loading drmaa support:\n{}".format(e))
        self.session = drmaa.Session()
        self.drmaa_args = drmaa_args
        self.drmaa_log_dir = drmaa_log_dir
        self.session.initialize()
        self.submitted = list()

    def cancel(self):
        from drmaa.const import JobControlAction
        from drmaa.errors import InvalidJobException, InternalException
        for jobid in self.submitted:
            try:
                self.session.control(jobid, JobControlAction.TERMINATE)
            except (InvalidJobException, InternalException):
                #This is common - logging a warning would probably confuse the user.
                pass
        self.shutdown()

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        super()._run(job)
        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        try:
            drmaa_args = job.format_wildcards(
                self.drmaa_args,
                cluster=self.cluster_wildcards(job))
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
        except (drmaa.InternalException,
                drmaa.InvalidAttributeValueException) as e:
            print_exception(WorkflowError("DRMAA Error: {}".format(e)),
                            self.workflow.linemaps)
            error_callback(job)
            return
        logger.info("Submitted DRMAA job {} with external jobid {}.".format(job.jobid, jobid))
        self.submitted.append(jobid)
        self.session.deleteJobTemplate(jt)

        submit_callback(job)

        with self.lock:
            self.active_jobs.append(DRMAAClusterJob(
                job, jobid, callback, error_callback, jobscript))

    def shutdown(self):
        super().shutdown()
        self.session.exit()

    def _wait_for_jobs(self):
        import drmaa
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
                        retval = self.session.wait(active_job.jobid,
                                                   drmaa.Session.TIMEOUT_NO_WAIT)
                    except drmaa.ExitTimeoutException as e:
                        # job still active
                        still_running.append(active_job)
                        continue
                    except (drmaa.InternalException, Exception) as e:
                        print_exception(WorkflowError("DRMAA Error: {}".format(e)),
                                        self.workflow.linemaps)
                        os.remove(active_job.jobscript)
                        active_job.error_callback(active_job.job)
                        continue
                    # job exited
                    os.remove(active_job.jobscript)
                    if retval.hasExited and retval.exitStatus == 0:
                        active_job.callback(active_job.job)
                    else:
                        self.print_job_error(active_job.job)
                        print_exception(
                            ClusterJobException(active_job, self.dag.jobid(active_job.job)),
                            self.workflow.linemaps)
                        active_job.error_callback(active_job.job)
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


@contextlib.contextmanager
def change_working_directory(directory=None):
    """ Change working directory in execution context if provided. """
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


KubernetesJob = namedtuple("KubernetesJob", "job jobid callback error_callback kubejob jobscript")


class KubernetesExecutor(ClusterExecutor):
    def __init__(self, workflow, dag, namespace, envvars,
                 container_image=None,
                 jobname="{rulename}.{jobid}",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 cluster_config=None,
                 local_input=None,
                 restart_times=None):

        exec_job = (
            'cp -rf /source/. . && '
            'snakemake {target} --snakefile {snakefile} '
            '--force -j{cores} --keep-target-files  --keep-remote '
            '--latency-wait 0 '
            ' --attempt {attempt} {use_threads} '
            '--wrapper-prefix {workflow.wrapper_prefix} '
            '{overwrite_config} {printshellcmds} --nocolor '
            '--notemp --no-hooks --nolock ')

        super().__init__(workflow, dag, None,
                         jobname=jobname,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         cluster_config=cluster_config,
                         local_input=local_input,
                         restart_times=restart_times,
                         exec_job=exec_job,
                         assume_shared_fs=False,
                         max_status_checks_per_second=10)
        # use relative path to Snakefile
        self.snakefile = os.path.relpath(workflow.snakefile)

        try:
            from kubernetes import config
        except ImportError:
            raise WorkflowError("The Python 3 package 'kubernetes' "
                                "must be installed to use Kubernetes")
        config.load_kube_config()

        import kubernetes.client
        self.kubeapi = kubernetes.client.CoreV1Api()
        self.batchapi = kubernetes.client.BatchV1Api()
        self.namespace = namespace
        self.envvars = envvars or []
        self.secret_files = {}
        self.run_namespace = str(uuid.uuid4())
        self.secret_envvars = {}
        self.register_secret()
        self.container_image = (
            container_image or
            get_container_image())

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
                logger.warning("Ignoring source file {}. Only files relative "
                               "to the working directory are allowed.".format(f))
                continue
            with open(f, "br") as content:
                key = "f{}".format(i)
                self.secret_files[key] = f
                secret.data[key] = base64.b64encode(content.read()).decode()
        for e in self.envvars:
            try:
                key = e.lower()
                secret.data[key] = base64.b64encode(os.environ[e].encode()).decode()
                self.secret_envvars[key] = e
            except KeyError:
                continue
        self.kubeapi.create_namespaced_secret(self.namespace, secret)

    def unregister_secret(self):
        import kubernetes.client
        self.kubeapi.delete_namespaced_secret(self.run_namespace,
                                              self.namespace,
                                              kubernetes.client.V1DeleteOptions())

    def shutdown(self):
        self.unregister_secret()
        super().shutdown()

    def cancel(self):
        import kubernetes.client
        body = kubernetes.client.V1DeleteOptions()
        with self.lock:
            for j in self.active_jobs:
                self.kubeapi.delete_namespaced_pod(
                    j.jobid, self.namespace, body)
        self.shutdown()

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        import kubernetes.client

        super()._run(job)
        exec_job = self.format_job(
            self.exec_job, job, _quote_all=True, rules=job.rules,
            use_threads="--force-use-threads" if not job.is_group() else "")
        # Kubernetes silently does not submit a job if the name is too long
        # therefore, we ensure that it is not longer than snakejob+uuid.
        jobid = "snakejob-{}".format(
            get_uuid("{}-{}-{}".format(
                self.run_namespace, job.jobid, job.attempt)))

        body = kubernetes.client.V1Pod()
        body.metadata = kubernetes.client.V1ObjectMeta(labels={"app": "snakemake"})
        body.metadata.name = jobid

        # container
        container = kubernetes.client.V1Container(name=jobid)
        container.image = self.container_image
        container.command = shlex.split("/bin/sh")
        container.args = ["-c", exec_job]
        container.working_dir = "/workdir"
        container.volume_mounts = [kubernetes.client.V1VolumeMount(
            name="workdir", mount_path="/workdir")]
        container.volume_mounts = [kubernetes.client.V1VolumeMount(
            name="source", mount_path="/source")]

        body.spec = kubernetes.client.V1PodSpec(containers=[container])
        # fail on first error
        body.spec.restart_policy = "Never"

        # source files as a secret volume
        # we copy these files to the workdir before executing Snakemake
        too_large = [path for path in self.secret_files.values()
                     if os.path.getsize(path) > 1000000]
        if too_large:
            raise WorkflowError("The following source files exceed the maximum "
                                "file size (1MB) that can be passed from host to "
                                "kubernetes. These are likely not source code "
                                "files. Consider adding them to your "
                                "remote storage instead or (if software) use "
                                "Conda packages or container images:\n{}".format(
                                "\n".join(too_large)))
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
                key=key, name=self.run_namespace)
            container.env.append(envvar)

        # request resources
        container.resources = kubernetes.client.V1ResourceRequirements()
        container.resources.requests = {}
        # Subtract 1 from the requested number of cores.
        # The reason is that kubernetes requires some cycles for
        # maintenance, but won't use a full core for that.
        # This way, we should be able to saturate the node without exceeding it
        # too much.
        container.resources.requests["cpu"] = job.resources["_cores"] - 1
        if "mem_mb" in job.resources.keys():
            container.resources.requests["memory"] = "{}M".format(
                job.resources["mem_mb"])

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
                privileged=True)

        pod = self._kubernetes_retry(
            lambda: self.kubeapi.create_namespaced_pod(self.namespace, body))

        logger.info("Get status with:\n"
                    "kubectl describe pod {jobid}\n"
                    "kubectl logs {jobid}".format(jobid=jobid))
        self.active_jobs.append(KubernetesJob(
            job, jobid, callback, error_callback, pod, None))

    def _kubernetes_retry(self, func):
        import kubernetes
        with self.lock:
            try:
                return func()
            except kubernetes.client.rest.ApiException as e:
                if e.status == 401:
                    # Unauthorized.
                    # Reload config in order to ensure token is
                    # refreshed. Then try again.
                    kubernetes.config.load_kube_config()
                    try:
                        return func()
                    except kubernetes.client.rest.ApiException as e:
                        # Both attempts failed, raise error.
                        raise WorkflowError(e,
                            "This is likely a bug in "
                            "https://github.com/kubernetes-client/python.")

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
                            lambda: self.kubeapi.read_namespaced_pod_status(j.jobid, self.namespace))
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
                        msg = ("Unknown pod {jobid}. "
                               "Has the pod been deleted "
                               "manually?").format(jobid=j.jobid)
                        self.print_job_error(j.job, msg=msg, jobid=j.jobid)
                        j.error_callback(j.job)
                    elif res.status.phase == "Failed":
                        msg = ("For details, please issue:\n"
                               "kubectl describe pod {jobid}\n"
                               "kubectl logs {jobid}").format(jobid=j.jobid)
                        # failed
                        self.print_job_error(j.job, msg=msg, jobid=j.jobid)
                        j.error_callback(j.job)
                    elif res.status.phase == "Succeeded":
                        # finished
                        j.callback(j.job)
                    else:
                        # still active
                        still_running.append(j)
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


def run_wrapper(job_rule, input, output, params, wildcards, threads, resources, log,
                benchmark, benchmark_repeats, conda_env, singularity_img,
                singularity_args, use_singularity, linemaps, debug,
                shadow_dir, jobid):
    """
    Wrapper around the run method that handles exceptions and benchmarking.

    Arguments
    job_rule   -- the ``job.rule`` member
    input      -- list of input files
    output     -- list of output files
    wildcards  -- so far processed wildcards
    threads    -- usable threads
    log        -- list of log files
    shadow_dir -- optional shadow directory root
    """
    # get shortcuts to job_rule members
    run = job_rule.run_func
    version = job_rule.version
    rule = job_rule.name
    is_shell = job_rule.shellcmd is not None

    if os.name == "posix" and debug:
        sys.stdin = open('/dev/stdin')

    if benchmark is not None:
        from snakemake.benchmark import BenchmarkRecord, benchmarked, write_benchmark_records

    try:
        with change_working_directory(shadow_dir):
            if benchmark:
                bench_records = []
                for bench_iteration in range(benchmark_repeats):
                    # Determine whether to benchmark this process or do not
                    # benchmarking at all.  We benchmark this process unless the
                    # execution is done through the ``shell:``, ``script:``, or
                    # ``wrapper:`` stanza.
                    is_sub = (job_rule.shellcmd or job_rule.script or
                              job_rule.wrapper or job_rule.cwl)
                    if is_sub:
                        # The benchmarking through ``benchmarked()`` is started
                        # in the execution of the shell fragment, script, wrapper
                        # etc, as the child PID is available there.
                        bench_record = BenchmarkRecord()
                        run(input, output, params, wildcards, threads, resources,
                            log, version, rule, conda_env, singularity_img,
                            singularity_args, use_singularity, bench_record,
                            jobid, is_shell, bench_iteration)
                    else:
                        # The benchmarking is started here as we have a run section
                        # and the generated Python function is executed in this
                        # process' thread.
                        with benchmarked() as bench_record:
                            run(input, output, params, wildcards, threads, resources,
                                log, version, rule, conda_env, singularity_img,
                                singularity_args, use_singularity,
                                bench_record, jobid, is_shell, bench_iteration)
                    # Store benchmark record for this iteration
                    bench_records.append(bench_record)
            else:
                run(input, output, params, wildcards, threads, resources,
                    log, version, rule, conda_env, singularity_img,
                    singularity_args, use_singularity, None, jobid, is_shell, None)
    except (KeyboardInterrupt, SystemExit) as e:
        # Re-raise the keyboard interrupt in order to record an error in the
        # scheduler but ignore it
        raise e
    except (Exception, BaseException) as ex:
        log_verbose_traceback(ex)
        # this ensures that exception can be re-raised in the parent thread
        lineno, file = get_exception_origin(ex, linemaps)
        raise RuleException(format_error(ex, lineno,
                                         linemaps=linemaps,
                                         snakefile=file,
                                         show_traceback=True))

    if benchmark is not None:
        try:
            write_benchmark_records(bench_records, benchmark)
        except (Exception, BaseException) as ex:
            raise WorkflowError(ex)
