__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
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
import tempfile
from functools import partial
from itertools import chain
from collections import namedtuple
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
from snakemake.common import Mode, __version__, get_container_image, get_uuid


# TODO move each executor into a separate submodule


def sleep():
    # do not sleep on CI. In that case we just want to quickly test everything.
    if os.environ.get("CI") != "true":
        time.sleep(10)


class AbstractExecutor:
    def __init__(
        self,
        workflow,
        dag,
        printreason=False,
        quiet=False,
        printshellcmds=False,
        printthreads=True,
        latency_wait=3,
        keepincomplete=False,
    ):
        self.workflow = workflow
        self.dag = dag
        self.quiet = quiet
        self.printreason = printreason
        self.printshellcmds = printshellcmds
        self.printthreads = printthreads
        self.latency_wait = latency_wait
        self.keepincomplete = keepincomplete

    def get_default_remote_provider_args(self):
        if self.workflow.default_remote_provider:
            return (
                " --default-remote-provider {} " "--default-remote-prefix {} "
            ).format(
                self.workflow.default_remote_provider.__module__.split(".")[-1],
                self.workflow.default_remote_prefix,
            )
        return ""

    def get_default_resources_args(self):
        if self.workflow.default_resources.args is not None:

            def fmt(res):
                if isinstance(res, str):
                    res = res.replace('"', r"\"")
                return '"{}"'.format(res)

            args = " --default-resources {} ".format(
                " ".join(map(fmt, self.workflow.default_resources.args))
            )
            return args
        return ""

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
        latency_wait=3,
        assume_shared_fs=True,
        keepincomplete=False,
    ):
        super().__init__(
            workflow,
            dag,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            keepincomplete=keepincomplete,
        )
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
        )
        self.stats.report_job_end(job)

    def handle_job_error(self, job, upload_remote=True):
        job.postprocess(
            error=True,
            assume_shared_fs=self.assume_shared_fs,
            latency_wait=self.latency_wait,
        )

    def get_additional_args(self):
        """Return a string to add to self.exec_job that includes additional
           arguments from the command line. This is currently used in the
           ClusterExecutor and CPUExecutor, as both were using the same
           code. Both have base class of the RealExecutor.
        """
        additional = ""
        if not self.workflow.cleanup_scripts:
            additional += " --skip-script-cleanup "
        if self.workflow.shadow_prefix:
            additional += " --shadow-prefix {} ".format(self.workflow.shadow_prefix)
        if self.workflow.use_conda:
            additional += " --use-conda "
            if self.workflow.conda_prefix:
                additional += " --conda-prefix {} ".format(self.workflow.conda_prefix)
        if self.workflow.use_singularity:
            additional += " --use-singularity "
            if self.workflow.singularity_prefix:
                additional += " --singularity-prefix {} ".format(
                    self.workflow.singularity_prefix
                )
            if self.workflow.singularity_args:
                additional += ' --singularity-args "{}"'.format(
                    self.workflow.singularity_args
                )

        if self.workflow.use_env_modules:
            additional += " --use-envmodules"

        return additional

    def format_job_pattern(self, pattern, job=None, **kwargs):
        overwrite_workdir = []
        if self.workflow.overwrite_workdir:
            overwrite_workdir.extend(("--directory", self.workflow.overwrite_workdir))

        overwrite_config = []
        if self.workflow.overwrite_configfiles:
            # add each of the overwriting configfiles in the original order
            if self.workflow.overwrite_configfiles:
                overwrite_config.append("--configfiles")
                overwrite_config.extend(self.workflow.overwrite_configfiles)
        if self.workflow.config_args:
            overwrite_config.append("--config")
            overwrite_config.extend(self.workflow.config_args)

        printshellcmds = ""
        if self.workflow.printshellcmds:
            printshellcmds = "-p"

        if not job.is_branched and not job.is_updated:
            # Restrict considered rules. This does not work for updated jobs
            # because they need to be updated in the spawned process as well.
            rules = ["--allowed-rules"]
            rules.extend(job.rules)
        else:
            rules = []

        target = kwargs.get("target", job.get_targets())
        snakefile = kwargs.get("snakefile", self.snakefile)
        cores = kwargs.get("cores", self.cores)
        if "target" in kwargs:
            del kwargs["target"]
        if "snakefile" in kwargs:
            del kwargs["snakefile"]
        if "cores" in kwargs:
            del kwargs["cores"]

        return format(
            pattern,
            job=job,
            attempt=job.attempt,
            overwrite_workdir=overwrite_workdir,
            overwrite_config=overwrite_config,
            printshellcmds=printshellcmds,
            workflow=self.workflow,
            snakefile=snakefile,
            cores=cores,
            benchmark_repeats=job.benchmark_repeats if not job.is_group() else None,
            target=target,
            rules=rules,
            **kwargs
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
        latency_wait=3,
        cores=1,
        keepincomplete=False,
    ):
        super().__init__(
            workflow,
            dag,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            keepincomplete=keepincomplete,
        )

        self.exec_job = "\\\n".join(
            (
                "cd {workflow.workdir_init} && ",
                "{sys.executable} -m snakemake {target} --snakefile {snakefile} ",
                "--force -j{cores} --keep-target-files --keep-remote ",
                "--attempt {attempt} ",
                "--force-use-threads --wrapper-prefix {workflow.wrapper_prefix} ",
                "--latency-wait {latency_wait} ",
                self.get_default_remote_provider_args(),
                self.get_default_resources_args(),
                "{overwrite_workdir} {overwrite_config} {printshellcmds} {rules} ",
                "--notemp --quiet --no-hooks --nolock --mode {} ".format(
                    Mode.subprocess
                ),
            )
        )

        self.exec_job += self.get_additional_args()
        self.use_threads = use_threads
        self.cores = cores

        # Zero thread jobs do not need a thread, but they occupy additional workers.
        # Hence we need to reserve additional workers for them.
        self.workers = workers + 5
        self.pool = concurrent.futures.ThreadPoolExecutor(max_workers=self.workers)

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

        conda_env = job.conda_env_path if self.workflow.use_conda else None
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
            self.workflow.edit_notebook,
        )

    def run_single_job(self, job):
        if self.use_threads or (not job.is_shadow and not job.is_run):
            future = self.pool.submit(
                self.cached_or_run, job, run_wrapper, *self.job_args_and_prepare(job)
            )
        else:
            # run directive jobs are spawned into subprocesses
            future = self.pool.submit(self.cached_or_run, job, self.spawn_job, job)
        return future

    def run_group_job(self, job):
        """Run a pipe group job.

        This lets all items run simultaneously."""
        # we only have to consider pipe groups because in local running mode,
        # these are the only groups that will occur

        futures = [self.run_single_job(j) for j in job]

        while True:
            k = 0
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
                        k += 1
            if k == len(futures):
                return
            time.sleep(1)

    def spawn_job(self, job):
        exec_job = self.exec_job
        cmd = self.format_job_pattern(
            exec_job, job=job, _quote_all=True, latency_wait=self.latency_wait
        )
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
        latency_wait=3,
        cluster_config=None,
        local_input=None,
        restart_times=None,
        exec_job=None,
        assume_shared_fs=True,
        max_status_checks_per_second=1,
        disable_default_remote_provider_args=False,
        disable_get_default_resources_args=False,
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
            latency_wait=latency_wait,
            assume_shared_fs=assume_shared_fs,
        )

        if not self.assume_shared_fs:
            # use relative path to Snakefile
            self.snakefile = os.path.relpath(workflow.snakefile)

        jobscript = workflow.jobscript
        if jobscript is None:
            jobscript = os.path.join(os.path.dirname(__file__), self.default_jobscript)
        try:
            with open(jobscript) as f:
                self.jobscript = f.read()
        except IOError as e:
            raise WorkflowError(e)

        if not "jobid" in get_wildcard_names(jobname):
            raise WorkflowError(
                'Defined jobname ("{}") has to contain the wildcard {jobid}.'
            )

        if exec_job is None:
            self.exec_job = "\\\n".join(
                (
                    "{envvars} " "cd {workflow.workdir_init} && "
                    if assume_shared_fs
                    else "",
                    "{sys.executable} " if assume_shared_fs else "python ",
                    "-m snakemake {target} --snakefile {snakefile} ",
                    "--force -j{cores} --keep-target-files --keep-remote ",
                    "--wait-for-files {wait_for_files} --latency-wait {latency_wait} ",
                    " --attempt {attempt} {use_threads} ",
                    "--wrapper-prefix {workflow.wrapper_prefix} ",
                    "{overwrite_workdir} {overwrite_config} {printshellcmds} {rules} "
                    "--nocolor --notemp --no-hooks --nolock ",
                    "--mode {} ".format(Mode.cluster),
                )
            )
        else:
            self.exec_job = exec_job

        self.exec_job += self.get_additional_args()
        if not disable_default_remote_provider_args:
            self.exec_job += self.get_default_remote_provider_args()
        if not disable_get_default_resources_args:
            self.exec_job += self.get_default_resources_args()

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
            max_calls=self.max_status_checks_per_second, period=1
        )

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

    def format_job(self, pattern, job, **kwargs):
        wait_for_files = []
        if self.assume_shared_fs:
            wait_for_files.append(self.tmpdir)
            wait_for_files.extend(job.get_wait_for_files())

        format_p = partial(
            self.format_job_pattern,
            job=job,
            properties=job.properties(cluster=self.cluster_params(job)),
            latency_wait=self.latency_wait,
            wait_for_files=wait_for_files,
            **kwargs
        )
        try:
            return format_p(pattern)
        except KeyError as e:
            raise WorkflowError(
                "Error formatting jobscript: {} not found\n"
                "Make sure that your custom jobscript is up to date.".format(e)
            )

    def write_jobscript(self, job, jobscript, **kwargs):
        # only force threads if this is not a group job
        # otherwise we want proper process handling
        use_threads = "--force-use-threads" if not job.is_group() else ""

        envvars = " ".join(
            "{}={}".format(var, os.environ[var]) for var in self.workflow.envvars
        )

        exec_job = self.format_job(
            self.exec_job,
            job,
            _quote_all=True,
            use_threads=use_threads,
            envvars=envvars,
            **kwargs
        )
        content = self.format_job(self.jobscript, job, exec_job=exec_job, **kwargs)
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
        cluster_config=None,
        jobname="snakejob.{rulename}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        latency_wait=3,
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
        self.external_jobid = dict()

        super().__init__(
            workflow,
            dag,
            cores,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            cluster_config=cluster_config,
            restart_times=restart_times,
            assume_shared_fs=assume_shared_fs,
            max_status_checks_per_second=max_status_checks_per_second,
        )

        if statuscmd:
            self.exec_job += " && exit 0 || exit 1"
        elif assume_shared_fs:
            # TODO wrap with watch and touch {jobrunning}
            # check modification date of {jobrunning} in the wait_for_job method
            self.exec_job += " && touch {jobfinished} || (touch {jobfailed}; exit 1)"
        else:
            raise WorkflowError(
                "If no shared filesystem is used, you have to "
                "specify a cluster status command."
            )

    def cancel(self):
        logger.info("Will exit after finishing currently running jobs.")
        self.shutdown()

    def register_job(self, job):
        # Do not register job here.
        # Instead do it manually once the jobid is known.
        pass

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        workdir = os.getcwd()
        jobid = job.jobid

        jobscript = self.get_jobscript(job)
        jobfinished = os.path.join(self.tmpdir, "{}.jobfinished".format(jobid))
        jobfailed = os.path.join(self.tmpdir, "{}.jobfailed".format(jobid))
        self.write_jobscript(
            job, jobscript, jobfinished=jobfinished, jobfailed=jobfailed
        )

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
            ext_jobid = (
                subprocess.check_output(
                    '{submitcmd} "{jobscript}"'.format(
                        submitcmd=submitcmd, jobscript=jobscript
                    ),
                    shell=True,
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
        if self.statuscmd is not None:

            def job_status(job):
                try:
                    # this command shall return "success", "failed" or "running"
                    return (
                        subprocess.check_output(
                            "{statuscmd} {jobid}".format(
                                jobid=job.jobid, statuscmd=self.statuscmd
                            ),
                            shell=True,
                        )
                        .decode()
                        .split("\n")[0]
                    )
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
                        raise WorkflowError(
                            "Failed to obtain job status. "
                            "See above for error message."
                        )

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
        latency_wait=3,
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
            latency_wait=latency_wait,
            cluster_config=cluster_config,
            restart_times=restart_times,
            assume_shared_fs=assume_shared_fs,
            max_status_checks_per_second=10,
        )
        self.submitcmd = submitcmd
        self.external_jobid = dict()

    def cancel(self):
        logger.info("Will exit after finishing currently running jobs.")
        self.shutdown()

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        workdir = os.getcwd()
        jobid = job.jobid

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
        latency_wait=3,
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
            latency_wait=latency_wait,
            cluster_config=cluster_config,
            restart_times=restart_times,
            assume_shared_fs=assume_shared_fs,
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
                        retval = self.session.wait(
                            active_job.jobid, drmaa.Session.TIMEOUT_NO_WAIT
                        )
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
                    # job exited
                    os.remove(active_job.jobscript)
                    if (
                        not retval.wasAborted
                        and retval.hasExited
                        and retval.exitStatus == 0
                    ):
                        active_job.callback(active_job.job)
                    else:
                        self.print_job_error(active_job.job)
                        self.print_cluster_job_error(
                            active_job, self.dag.jobid(active_job.job)
                        )
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
        latency_wait=3,
        cluster_config=None,
        local_input=None,
        restart_times=None,
        keepincomplete=False,
    ):

        exec_job = (
            "cp -rf /source/. . && "
            "snakemake {target} --snakefile {snakefile} "
            "--force -j{cores} --keep-target-files  --keep-remote "
            "--latency-wait 0 "
            " --attempt {attempt} {use_threads} "
            "--wrapper-prefix {workflow.wrapper_prefix} "
            "{overwrite_config} {printshellcmds} {rules} --nocolor "
            "--notemp --no-hooks --nolock "
        )

        super().__init__(
            workflow,
            dag,
            None,
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
            max_status_checks_per_second=10,
        )
        # use relative path to Snakefile
        self.snakefile = os.path.relpath(workflow.snakefile)

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
            max_file_size = 1000000
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

        self.kubeapi.delete_namespaced_secret(
            self.run_namespace, self.namespace, body=kubernetes.client.V1DeleteOptions()
        )

    def shutdown(self):
        self.unregister_secret()
        super().shutdown()

    def cancel(self):
        import kubernetes.client

        body = kubernetes.client.V1DeleteOptions()
        with self.lock:
            for j in self.active_jobs:
                self.kubeapi.delete_namespaced_pod(j.jobid, self.namespace, body=body)
        self.shutdown()

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        import kubernetes.client

        super()._run(job)
        exec_job = self.format_job(
            self.exec_job,
            job,
            _quote_all=True,
            use_threads="--force-use-threads" if not job.is_group() else "",
        )
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

        body.spec = kubernetes.client.V1PodSpec(containers=[container])
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
                    logger.info("trying to reauthenticate")
                    kubernetes.config.load_kube_config()
                    subprocess.run(["kubectl", "get", "nodes"])

                    self.kubeapi = kubernetes.client.CoreV1Api()
                    self.batchapi = kubernetes.client.BatchV1Api()
                    self.register_secret()
                    try:
                        return func()
                    except kubernetes.client.rest.ApiException as e:
                        # Both attempts failed, raise error.
                        raise WorkflowError(
                            e,
                            "This is likely a bug in "
                            "https://github.com/kubernetes-client/python.",
                        )
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
                        body = kubernetes.client.V1DeleteOptions()
                        self.kubeapi.delete_namespaced_pod(
                            j.jobid, self.namespace, body=body
                        )
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
        latency_wait=3,
        local_input=None,
        restart_times=None,
        max_status_checks_per_second=1,
        keepincomplete=False,
    ):

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
        self.snakefile = workflow.snakefile
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
        exec_job = (
            "snakemake {target} --snakefile {snakefile} "
            "--force -j{cores} --keep-target-files  --keep-remote "
            "--latency-wait 0 "
            "--attempt 1 {use_threads} "
            "{overwrite_config} {rules} --nocolor "
            "--notemp --no-hooks --nolock "
        )

        super().__init__(
            workflow,
            dag,
            cores,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            local_input=local_input,
            restart_times=restart_times,
            exec_job=exec_job,
            assume_shared_fs=False,
            max_status_checks_per_second=max_status_checks_per_second,
            disable_default_remote_provider_args=True,
            disable_get_default_resources_args=True,
        )
        self.container_image = container_image or get_container_image()
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

    def handle_remote(self, target):
        if isinstance(target, _IOFile) and target.remote_object.provider.is_default:
            return self.remove_prefix(target)
        else:
            return target

    def add_command(self, job, tibanna_args, tibanna_config):
        # snakefile, with file name remapped
        snakefile_fname = tibanna_args.snakemake_main_filename
        # targets, with file name remapped
        targets = job.get_targets()
        if not isinstance(targets, list):
            targets = [targets]
        targets_default = " ".join([self.handle_remote(t) for t in targets])
        # use_threads
        use_threads = "--force-use-threads" if not job.is_group() else ""
        # format command
        command = self.format_job_pattern(
            self.exec_job,
            job,
            target=targets_default,
            snakefile=snakefile_fname,
            use_threads=use_threads,
            cores=tibanna_config["cpu"],
        )
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
):
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
                )
    except (KeyboardInterrupt, SystemExit) as e:
        # Re-raise the keyboard interrupt in order to record an error in the
        # scheduler but ignore it
        raise e
    except (Exception, BaseException) as ex:
        log_verbose_traceback(ex)
        # this ensures that exception can be re-raised in the parent thread
        lineno, file = get_exception_origin(ex, linemaps)
        raise RuleException(
            format_error(
                ex, lineno, linemaps=linemaps, snakefile=file, show_traceback=True
            )
        )

    if benchmark is not None:
        try:
            write_benchmark_records(bench_records, benchmark)
        except (Exception, BaseException) as ex:
            raise WorkflowError(ex)
