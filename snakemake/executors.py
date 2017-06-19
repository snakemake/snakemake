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
import threading
import concurrent.futures
import subprocess
import signal
from functools import partial
from itertools import chain
from collections import namedtuple
from tempfile import mkdtemp

from snakemake.jobs import Job
from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.stats import Stats
from snakemake.utils import format, Unformattable, makedirs
from snakemake.io import get_wildcard_names, Wildcards
from snakemake.exceptions import print_exception, get_exception_origin
from snakemake.exceptions import format_error, RuleException, log_verbose_traceback
from snakemake.exceptions import ClusterJobException, ProtectedOutputException, WorkflowError, ImproperShadowException, SpawnedJobError
from snakemake.common import Mode


class AbstractExecutor:
    def __init__(self, workflow, dag,
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 printthreads=True,
                 latency_wait=3,
                 benchmark_repeats=1):
        self.workflow = workflow
        self.dag = dag
        self.quiet = quiet
        self.printreason = printreason
        self.printshellcmds = printshellcmds
        self.printthreads = printthreads
        self.latency_wait = latency_wait
        self.benchmark_repeats = benchmark_repeats

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        job.check_protected_output()
        self._run(job)
        callback(job)

    def shutdown(self):
        pass

    def cancel(self):
        pass

    def _run(self, job):
        self.printjob(job)

    def rule_prefix(self, job):
        return "local " if self.workflow.is_local(job.rule) else ""

    def printjob(self, job):
        # skip dynamic jobs that will be "executed" only in dryrun mode
        if self.dag.dynamic(job):
            return

        def format_files(job, io, dynamicio):
            for f in io:
                if f in dynamicio:
                    yield "{} (dynamic)".format(f.format_dynamic())
                else:
                    yield f

        priority = self.dag.priority(job)
        logger.job_info(jobid=self.dag.jobid(job),
                        msg=job.message,
                        name=job.rule.name,
                        local=self.workflow.is_local(job.rule),
                        input=list(format_files(job, job.input,
                                                job.dynamic_input)),
                        output=list(format_files(job, job.output,
                                                 job.dynamic_output)),
                        log=list(job.log),
                        benchmark=job.benchmark,
                        wildcards=job.wildcards_dict,
                        reason=str(self.dag.reason(job)),
                        resources=job.resources,
                        priority="highest"
                        if priority == Job.HIGHEST_PRIORITY else priority,
                        threads=job.threads)

        if job.dynamic_output:
            logger.info("Subsequent jobs will be added dynamically "
                        "depending on the output of this rule")

    def print_job_error(self, job):
        logger.error("Error in job {} while creating output file{} {}.".format(
            job, "s" if len(job.output) > 1 else "", ", ".join(job.output)))

    def handle_job_success(self, job):
        pass

    def handle_job_error(self, job):
        pass


class DryrunExecutor(AbstractExecutor):
    def _run(self, job):
        super()._run(job)
        logger.shellcmd(job.shellcmd)


class RealExecutor(AbstractExecutor):
    def __init__(self, workflow, dag,
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 benchmark_repeats=1):
        super().__init__(workflow, dag,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         benchmark_repeats=benchmark_repeats)
        self.stats = Stats()

    def _run(self, job, callback=None, error_callback=None):
        super()._run(job)
        self.stats.report_job_start(job)
        try:
            self.workflow.persistence.started(job)
        except IOError as e:
            logger.info(
                "Failed to set marker file for job started ({}). "
                "Snakemake will work, but cannot ensure that output files "
                "are complete in case of a kill signal or power loss. "
                "Please ensure write permissions for the "
                "directory {}".format(e, self.workflow.persistence.path))

    def handle_job_success(self, job, upload_remote=True, ignore_missing_output=False):
        self.dag.handle_touch(job)
        self.dag.check_and_touch_output(
            job,
            wait=self.latency_wait,
            ignore_missing_output=ignore_missing_output)
        self.dag.unshadow_output(job)
        self.dag.handle_remote(job, upload=upload_remote)
        self.dag.handle_protected(job)
        self.dag.handle_temp(job)
        job.close_remote()

        self.stats.report_job_end(job)
        try:
            self.workflow.persistence.finished(job)
        except IOError as e:
            logger.info("Failed to remove marker file for job started "
                        "({}). Please ensure write permissions for the "
                        "directory {}".format(e,
                                              self.workflow.persistence.path))

    def handle_job_error(self, job):
        job.close_remote()

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

        target = job.output if job.output else job.rule.name

        return format(pattern,
                      job=job,
                      overwrite_workdir=overwrite_workdir,
                      overwrite_config=overwrite_config,
                      printshellcmds=printshellcmds,
                      workflow=self.workflow,
                      cores=self.cores,
                      benchmark_repeats=self.benchmark_repeats,
                      target=target,
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
                 benchmark_repeats=1,
                 cores=1):
        super().__init__(workflow, dag,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         benchmark_repeats=benchmark_repeats)

        self.exec_job = '\\\n'.join((
            'cd {workflow.workdir_init} && ',
            '{sys.executable} -m snakemake {target} --snakefile {workflow.snakefile} ',
            '--force -j{cores} --keep-target-files --keep-shadow --keep-remote ',
            '--benchmark-repeats {benchmark_repeats} ',
            '--force-use-threads --wrapper-prefix {workflow.wrapper_prefix} ',
            '{overwrite_workdir} {overwrite_config} {printshellcmds} ',
            '--notemp --quiet --no-hooks --nolock --mode {} '.format(Mode.subprocess)))

        if self.workflow.use_conda:
            self.exec_job += " --use-conda "
            if self.workflow.conda_prefix:
                self.exec_job += " --conda-prefix " + self.workflow.conda_prefix + " "

        self.use_threads = use_threads
        self.cores = cores
        self.pool = concurrent.futures.ThreadPoolExecutor(max_workers=workers)

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        super()._run(job)

        if self.use_threads or (not job.is_shadow and (job.is_shell or job.is_norun or job.is_script or job.is_wrapper)):
            job.prepare()
            conda_env = None
            if self.workflow.use_conda:
                conda_env = job.conda_env

            benchmark = None
            if job.benchmark is not None:
                benchmark = str(job.benchmark)
            future = self.pool.submit(
                run_wrapper, job.rule, job.input.plainstrings(),
                job.output.plainstrings(), job.params, job.wildcards, job.threads,
                job.resources, job.log.plainstrings(), benchmark,
                self.benchmark_repeats, conda_env,
                self.workflow.linemaps, self.workflow.debug,
                shadow_dir=job.shadow_dir)
        else:
            # run directive jobs are spawned into subprocesses
            future = self.pool.submit(self.spawn_job, job)

        future.add_done_callback(partial(self._callback, job, callback,
                                         error_callback))

    def spawn_job(self, job):
        exec_job = self.exec_job
        if not job.rule.is_branched:
            exec_job += " --allowed-rules {}".format(job.rule)
        cmd = self.format_job_pattern(exec_job, job=job, _quote_all=True)
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
            if ex:
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
        job.cleanup()
        self.workflow.persistence.cleanup(job)


class ClusterExecutor(RealExecutor):

    default_jobscript = "jobscript.sh"

    def __init__(self, workflow, dag, cores,
                 jobname="snakejob.{rulename}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 benchmark_repeats=1,
                 cluster_config=None,
                 local_input=None,
                 max_jobs_per_second=None,
                 restart_times=None):
        local_input = local_input or []
        super().__init__(workflow, dag,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         benchmark_repeats=benchmark_repeats)

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

        self.exec_job = '\\\n'.join((
            'cd {workflow.workdir_init} && ',
            '{sys.executable} -m snakemake {target} --snakefile {workflow.snakefile} ',
            '--force -j{cores} --keep-target-files --keep-shadow --keep-remote ',
            '--wait-for-files {wait_for_files} --latency-wait {latency_wait} ',
            '--benchmark-repeats {benchmark_repeats} ',
            '--force-use-threads --wrapper-prefix {workflow.wrapper_prefix} ',
            '{overwrite_workdir} {overwrite_config} {printshellcmds} --nocolor ',
            '--notemp --quiet --no-hooks --nolock'))

        if printshellcmds:
            self.exec_job += " --printshellcmds "
        if self.workflow.use_conda:
            self.exec_job += " --use-conda "
            if self.workflow.conda_prefix:
                self.exec_job += " --conda-prefix " + self.workflow.conda_prefix + " "

        # force threading.Lock() for cluster jobs
        self.exec_job += " --force-use-threads "

        if not any(dag.dynamic_output_jobs):
            # disable restiction to target rule in case of dynamic rules!
            self.exec_job += " --allowed-rules {job.rule.name} "
        self.jobname = jobname
        self._tmpdir = None
        self.cores = cores if cores else ""
        self.cluster_config = cluster_config if cluster_config else dict()

        self.max_jobs_per_second = max_jobs_per_second
        self.restart_times = restart_times
        if self.max_jobs_per_second:
            self.rate_lock = threading.RLock()
            self.rate_interval = 1 / self.max_jobs_per_second
            self.rate_last_called = 0

        self.active_jobs = list()
        self.lock = threading.Lock()
        self.wait = True
        self.wait_thread = threading.Thread(target=self._wait_for_jobs)
        self.wait_thread.daemon = True
        self.wait_thread.start()

    def shutdown(self):
        with self.lock:
            self.wait = False
        self.wait_thread.join()
        shutil.rmtree(self.tmpdir)

    def cancel(self):
        self.shutdown()

    def _limit_rate(self):
        """Called in ``_run()`` for rate-limiting"""
        with self.rate_lock:
            elapsed = time.clock() - self.rate_last_called
            wait = self.rate_interval - elapsed
            if wait > 0:
                time.sleep(wait)
            self.rate_last_called = time.clock()

    def _run(self, job, callback=None, error_callback=None):
        if self.max_jobs_per_second:
            self._limit_rate()
        job.remove_existing_output()
        job.download_remote_input()
        super()._run(job, callback=callback, error_callback=error_callback)
        logger.shellcmd(job.shellcmd)

    @property
    def tmpdir(self):
        if self._tmpdir is None:
            self._tmpdir = mkdtemp(dir=".snakemake", prefix="tmp.")
        return os.path.abspath(self._tmpdir)

    def get_jobscript(self, job):
        f = job.format_wildcards(self.jobname,
                             rulename=job.rule.name,
                             jobid=self.dag.jobid(job),
                             cluster=self.cluster_wildcards(job))
        if os.path.sep in f:
            raise WorkflowError("Path separator ({}) found in job name {}. "
                                "This is not supported.".format(
                                os.path.sep, f))

        return os.path.join(self.tmpdir, f)

    def spawn_jobscript(self, job, jobscript, **kwargs):
        wait_for_files = [self.tmpdir]
        wait_for_files.extend(job.local_input)
        wait_for_files.extend(f.local_file()
                              for f in job.remote_input if not f.stay_on_remote)

        if job.shadow_dir:
            wait_for_files.append(job.shadow_dir)
        if self.workflow.use_conda and job.conda_env:
            wait_for_files.append(job.conda_env)

        format_p = partial(self.format_job_pattern,
                           job=job,
                           properties=json.dumps(job.properties(
                               cluster=self.cluster_params(job))),
                           latency_wait=self.latency_wait,
                           wait_for_files=wait_for_files,
                           **kwargs)
        exec_job = self.exec_job
        try:
            exec_job = format_p(exec_job, _quote_all=True)
            with open(jobscript, "w") as f:
                print(format_p(self.jobscript, exec_job=exec_job), file=f)
        except KeyError as e:
            raise WorkflowError(
                "Error formatting jobscript: {} not found\n"
                "Make sure that your custom jobscript is up to date.".format(e))
        os.chmod(jobscript, os.stat(jobscript).st_mode | stat.S_IXUSR)

    def cluster_params(self, job):
        """Return wildcards object for job from cluster_config."""

        cluster = self.cluster_config.get("__default__", dict()).copy()
        cluster.update(self.cluster_config.get(job.rule.name, dict()))
        # Format values with available parameters from the job.
        for key, value in list(cluster.items()):
            if isinstance(value, str):
                cluster[key] = job.format_wildcards(value)

        return cluster

    def cluster_wildcards(self, job):
        return Wildcards(fromdict=self.cluster_params(job))

    def handle_job_success(self, job):
        super().handle_job_success(job, upload_remote=False)

    def handle_job_error(self, job):
        # TODO what about removing empty remote dirs?? This cannot be decided
        # on the cluster node.
        super().handle_job_error(job)


GenericClusterJob = namedtuple("GenericClusterJob", "job jobid callback error_callback jobscript jobfinished jobfailed")


class GenericClusterExecutor(ClusterExecutor):
    def __init__(self, workflow, dag, cores,
                 submitcmd="qsub",
                 cluster_config=None,
                 jobname="snakejob.{rulename}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 benchmark_repeats=1,
                 max_jobs_per_second=None,
                 restart_times=0):
        super().__init__(workflow, dag, cores,
                         jobname=jobname,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         benchmark_repeats=benchmark_repeats,
                         cluster_config=cluster_config,
                         max_jobs_per_second=max_jobs_per_second,
                         restart_times=restart_times)
        self.submitcmd = submitcmd
        self.external_jobid = dict()
        # TODO wrap with watch and touch {jobrunning}
        # check modification date of {jobrunning} in the wait_for_job method
        self.exec_job += ' && touch "{jobfinished}" || (touch "{jobfailed}"; exit 1)'

    def cancel(self):
        logger.info("Will exit after finishing currently running jobs.")
        self.shutdown()

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        super()._run(job)
        workdir = os.getcwd()
        jobid = self.dag.jobid(job)

        jobscript = self.get_jobscript(job)
        jobfinished = os.path.join(self.tmpdir, "{}.jobfinished".format(jobid))
        jobfailed = os.path.join(self.tmpdir, "{}.jobfailed".format(jobid))
        self.spawn_jobscript(job, jobscript,
                             jobfinished=jobfinished,
                             jobfailed=jobfailed)

        deps = " ".join(self.external_jobid[f] for f in job.input
                        if f in self.external_jobid)
        try:
            submitcmd = job.format_wildcards(
                self.submitcmd,
                dependencies=deps,
                cluster=self.cluster_wildcards(job))
        except AttributeError as e:
            raise WorkflowError(str(e), rule=job.rule)
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
            logger.debug("Submitted job {} with external jobid {}.".format(
                jobid, ext_jobid))

        submit_callback(job)
        with self.lock:
            self.active_jobs.append(GenericClusterJob(job, ext_jobid, callback, error_callback, jobscript, jobfinished, jobfailed))

    def _wait_for_jobs(self):
        while True:
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                for active_job in active_jobs:
                    if os.path.exists(active_job.jobfinished):
                        os.remove(active_job.jobfinished)
                        os.remove(active_job.jobscript)
                        active_job.callback(active_job.job)
                    elif os.path.exists(active_job.jobfailed):
                        os.remove(active_job.jobfailed)
                        os.remove(active_job.jobscript)
                        self.print_job_error(active_job.job)
                        print_exception(ClusterJobException(active_job, self.dag.jobid(active_job.job)),
                                        self.workflow.linemaps)
                        active_job.error_callback(active_job.job)
                    else:
                        self.active_jobs.append(active_job)
            time.sleep(1)


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
                 benchmark_repeats=1,
                 max_jobs_per_second=None,
                 restart_times=0):
        super().__init__(workflow, dag, cores,
                         jobname=jobname,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         benchmark_repeats=benchmark_repeats,
                         cluster_config=cluster_config,
                         max_jobs_per_second=max_jobs_per_second,
                         restart_times=restart_times)
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
        jobid = self.dag.jobid(job)

        jobscript = self.get_jobscript(job)
        self.spawn_jobscript(job, jobscript)

        deps = " ".join(self.external_jobid[f] for f in job.input
                        if f in self.external_jobid)
        try:
            submitcmd = job.format_wildcards(
                self.submitcmd,
                dependencies=deps,
                cluster=self.cluster_wildcards(job))
        except AttributeError as e:
            raise WorkflowError(str(e), rule=job.rule)

        process = subprocess.Popen('{submitcmd} "{jobscript}"'.format(submitcmd=submitcmd,
                                           jobscript=jobscript), shell=True)
        submit_callback(job)

        with self.lock:
            self.active_jobs.append(SynchronousClusterJob(job, process.pid, callback, error_callback, jobscript, process))

    def _wait_for_jobs(self):
        while True:
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                for active_job in active_jobs:
                    exitcode = active_job.process.poll()
                    if exitcode is None:
                        # job not yet finished
                        self.active_jobs.append(active_job)
                    elif exitcode == 0:
                        # job finished successfully
                        os.remove(active_job.jobscript)
                        active_job.callback(active_job.job)
                    else:
                        # job failed
                        os.remove(active_job.jobscript)
                        self.print_job_error(active_job.job)
                        print_exception(ClusterJobException(active_job, self.dag.jobid(active_job.job)),
                                        self.workflow.linemaps)
                        active_job.error_callback(active_job.job)
            time.sleep(1)


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
                 benchmark_repeats=1,
                 cluster_config=None,
                 max_jobs_per_second=None,
                 restart_times=0):
        super().__init__(workflow, dag, cores,
                         jobname=jobname,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         latency_wait=latency_wait,
                         benchmark_repeats=benchmark_repeats,
                         cluster_config=cluster_config,
                         max_jobs_per_second=max_jobs_per_second,
                         restart_times=restart_times)
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
        self.spawn_jobscript(job, jobscript)

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
        logger.info("Submitted DRMAA job {} with external jobid {}.".format(self.dag.jobid(job), jobid))
        self.submitted.append(jobid)
        self.session.deleteJobTemplate(jt)

        submit_callback(job)

        with self.lock:
            self.active_jobs.append(DRMAAClusterJob(job, jobid, callback, error_callback, jobscript))

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
                for active_job in active_jobs:
                    try:
                        retval = self.session.wait(active_job.jobid,
                                                   drmaa.Session.TIMEOUT_NO_WAIT)
                    except drmaa.ExitTimeoutException as e:
                        # job still active
                        self.active_jobs.append(active_job)
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
            time.sleep(1)


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


def run_wrapper(job_rule, input, output, params, wildcards, threads, resources, log,
                benchmark, benchmark_repeats, conda_env, linemaps, debug=False,
                shadow_dir=None):
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

    if os.name == "posix" and debug:
        sys.stdin = open('/dev/stdin')

    if benchmark is not None:
        from snakemake.benchmark import BenchmarkRecord, benchmarked, write_benchmark_records

    try:
        with change_working_directory(shadow_dir):
            if benchmark:
                bench_records = []
                for i in range(benchmark_repeats):
                    # Determine whether to benchmark this process or do not
                    # benchmarking at all.  We benchmark this process unless the
                    # execution is done through the ``shell:``, ``script:``, or
                    # ``wrapper:`` stanza.
                    is_sub = job_rule.shellcmd or job_rule.script or job_rule.wrapper
                    if is_sub:
                        # The benchmarking through ``benchmarked()`` is started
                        # in the execution of the shell fragment, script, wrapper
                        # etc, as the child PID is available there.
                        bench_record = BenchmarkRecord()
                        run(input, output, params, wildcards, threads, resources,
                            log, version, rule, conda_env, bench_record)
                    else:
                        # The benchmarking is started here as we have a run section
                        # and the generated Python function is executed in this
                        # process' thread.
                        with benchmarked() as bench_record:
                            run(input, output, params, wildcards, threads, resources,
                                log, version, rule, conda_env, bench_record)
                    # Store benchmark record for this iteration
                    bench_records.append(bench_record)
            else:
                run(input, output, params, wildcards, threads, resources,
                    log, version, rule, conda_env, None)
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
