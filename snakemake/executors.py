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
import tarfile
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
from snakemake.common import (
    bytesto,
    Mode,
    __version__,
    get_container_image,
    get_uuid,
    get_file_hash,
)


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
            args = " --default-resources {} ".format(
                " ".join(map('"{}"'.format, self.workflow.default_resources.args))
            )
            return args
        return ""

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
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


GoogleLifeSciencesJob = namedtuple(
    "GoogleLifeSciencesJob", "job jobname jobid callback error_callback"
)


class GoogleLifeSciencesExecutor(ClusterExecutor):
    """the GoogleLifeSciences executor uses Google Cloud Storage, and
       Compute Engine paired with the Google Life Sciences API.
       https://cloud.google.com/life-sciences/docs/quickstart
    """

    def __init__(
        self,
        workflow,
        dag,
        cores,
        envvars,
        jobname="snakejob.{name}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        container_image=None,
        regions=None,
        location=None,
        cache=False,
        machine_type_prefix=None,
        latency_wait=3,
        local_input=None,
        restart_times=None,
        exec_job=None,
        max_status_checks_per_second=1,
    ):

        # Attach variables for easy access
        self.workflow = workflow
        self.quiet = quiet
        self.workdir = os.path.dirname(self.workflow.persistence.path)
        self._save_storage_cache = cache
        self._machine_type_prefix = machine_type_prefix

        # Relative path for running on instance
        self._set_snakefile()

        # Prepare workflow sources for build package
        self._set_workflow_sources()

        exec_job = exec_job or (
            "snakemake {target} --snakefile %s "
            "--force -j{cores} --keep-target-files --keep-remote "
            "--latency-wait 0 "
            "--attempt 1 {use_threads} "
            "{overwrite_config} {rules} --nocolor "
            "--notemp --no-hooks --nolock " % self.snakefile
        )

        # IMPORTANT: using Compute Engine API and not k8s == no support secrets
        self.envvars = envvars or []

        # Quit early if we can't authenticate
        self._get_services()
        self._get_bucket()

        # Akin to Kubernetes, create a run namespace, default container image
        self.run_namespace = str(uuid.uuid4())
        self.container_image = container_image or get_container_image()
        self.regions = regions or ["us-east1", "us-west1", "us-central1"]

        # Determine API location based on user preference, and then regions
        self._set_location(location)

        # Tell the user right away the regions, location, and container
        logger.debug("regions=%s" % self.regions)
        logger.debug("location=%s" % self.location)
        logger.debug("container=%s" % self.container_image)

        # The project name is required, either from client or environment
        self.project = (
            os.environ.get("GOOGLE_CLOUD_PROJECT") or self._bucket_service.project
        )
        if not self.project:
            raise WorkflowError(
                "You must provide a --google-project or export "
                "GOOGLE_CLOUD_PROJECT to use the Life Sciences API. "
            )

        # Keep track of build packages to clean up shutdown
        self._build_packages = set()

        # Save default resources to add later, since we need to add custom
        # default resources depending on the instance requested
        self.default_resources = self.workflow.default_resources
        self.workflow.default_resources.args = None

        super().__init__(
            workflow,
            dag,
            None,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            restart_times=restart_times,
            exec_job=exec_job,
            assume_shared_fs=False,
            max_status_checks_per_second=10,
        )

    def _get_services(self):
        """use the Google Discovery Build to generate API clients
           for Life Sciences, and use the google storage python client 
           for storage.
        """
        from googleapiclient.discovery import build as discovery_build
        from oauth2client.client import (
            GoogleCredentials,
            ApplicationDefaultCredentialsError,
        )
        from google.cloud import storage

        # Credentials must be exported to environment
        try:
            creds = GoogleCredentials.get_application_default()
        except ApplicationDefaultCredentialsError as ex:
            log_verbose_traceback(ex)
            raise ex

        # Discovery clients for Google Cloud Storage and Life Sciences API
        self._storage_cli = discovery_build("storage", "v1", credentials=creds)
        self._compute_cli = discovery_build("compute", "v1", credentials=creds)
        self._api = discovery_build("lifesciences", "v2beta", credentials=creds)
        self._bucket_service = storage.Client()

    def _get_bucket(self):
        """get a connection to the storage bucket (self.bucket) and exit
           if the name is taken or otherwise invalid.

           Parameters
           ==========
           workflow: the workflow object to derive the prefix from
        """
        import google

        # Hold path to requested subdirectory and main bucket
        bucket_name = self.workflow.default_remote_prefix.split("/")[0]
        self.gs_subdir = re.sub(
            "^{}/".format(bucket_name), "", self.workflow.default_remote_prefix
        )

        # Case 1: The bucket already exists
        try:
            self.bucket = self._bucket_service.get_bucket(bucket_name)

        # Case 2: The bucket needs to be created
        except google.cloud.exceptions.NotFound:
            self.bucket = self._bucket_service.create_bucket(bucket_name)

        # Case 2: The bucket name is already taken
        except (Exception, BaseException) as ex:
            logger.error(
                "Cannot get or create {} (exit code {}):\n{}".format(
                    bucket_name, ex.returncode, ex.output.decode()
                )
            )
            log_verbose_traceback(ex)
            raise ex

        logger.debug("bucket=%s" % self.bucket.name)
        logger.debug("subdir=%s" % self.gs_subdir)

    def _set_location(self, location=None):
        """The location is where the Google Life Sciences API is located.
           This can be meaningful if the requester has data residency
           requirements or multi-zone needs. To determine this value,
           we first use the locations API to determine locations available,
           and then compare them against:

           1. user specified location or prefix
           2. regions having the same prefix
           3. if cannot be satisifed, we throw an error.
        """
        # Derive available locations
        # See https://cloud.google.com/life-sciences/docs/concepts/locations
        locations = (
            self._api.projects()
            .locations()
            .list(name="projects/snakemake-testing")
            .execute()
        )

        locations = {x["locationId"]: x["name"] for x in locations.get("locations", [])}

        # Alert the user about locations available
        logger.debug("locations-available:\n%s" % "\n".join(locations))

        # If no locations, there is something wrong
        if not locations:
            raise WorkflowError("No locations found for Google Life Sciences API.")

        # First pass, attempt to match the user-specified location (or prefix)
        if location:
            if location in locations:
                self.location = locations[location]
                return

            # It could be that a prefix was provided
            for contender in locations:
                if contender.startswith(location):
                    self.location = locations[contender]
                    return

            # If we get here and no match, alert user.
            raise WorkflowError(
                "Location or prefix requested %s is not available." % location
            )

        # If we get here, we need to select location from regions
        for region in self.regions:
            if region in locations:
                self.location = locations[region]
                return

        # If we get here, choose based on prefix
        prefixes = set([r.split("-")[0] for r in self.regions])
        regexp = "^(%s)" % "|".join(prefixes)
        for location in locations:
            if re.search(regexp, location):
                self.location = locations[location]
                return

        # If we get here, total failure of finding location
        raise WorkflowError(
            " No locations available for regions!"
            " Please specify a location with --google-lifesciences-location "
            " or extend --google-lifesciences-regions to find a Life Sciences location."
        )

    def shutdown(self):
        """shutdown deletes build packages if the user didn't request to clean
           up the cache. At this point we've already cancelled running jobs.
        """
        # Delete build source packages only if user regooglquested no cache
        if self._save_storage_cache:
            logger.debug("Requested to save workflow sources, skipping cleanup.")
        else:
            for package in self._build_packages:
                blob = self.bucket.blob(package)
                if blob.exists():
                    logger.debug("Deleting blob %s" % package)
                    blob.delete()

        # perform additional steps on shutdown if necessary
        super().shutdown()

    def cancel(self):
        """cancel execution, usually by way of control+c. Cleanup is done in
           shutdown (deleting cached workdirs in Google Cloud Storage
        """
        import googleapiclient

        # projects.locations.operations/cancel
        operations = self._api.projects().locations().operations()

        for job in self.active_jobs:
            request = operations.cancel(name=job.jobname)
            logger.debug("Cancelling operation {}".format(job.jobid))
            try:
                self._retry_request(request)
            except (Exception, BaseException, googleapiclient.errors.HttpError):
                continue

        self.shutdown()

    def get_available_machine_types(self):
        """Using the regions available at self.regions, use the GCP API
           to retrieve a lookup dictionary of all available machine types.
        """
        # Regular expression to determine if zone in region
        regexp = "^(%s)" % "|".join(self.regions)

        # Retrieve zones, filter down to selected regions
        zones = self._retry_request(
            self._compute_cli.zones().list(project=self.project)
        )
        zones = [z for z in zones["items"] if re.search(regexp, z["name"])]

        # Retrieve machine types available across zones
        # https://cloud.google.com/compute/docs/regions-zones/
        lookup = {}
        for zone in zones:
            request = self._compute_cli.machineTypes().list(
                project=self.project, zone=zone["name"]
            )
            lookup[zone["name"]] = self._retry_request(request)["items"]

        # Only keep those that are shared, use last zone as a base
        machine_types = {mt["name"]: mt for mt in lookup[zone["name"]]}
        del lookup[zone["name"]]

        # Update final list based on the remaining
        to_remove = set()
        for zone, types in lookup.items():
            names = [x["name"] for x in types]
            names = [name for name in names if "micro" not in name]
            names = [name for name in names if not re.search("^(e2|m1)", name)]
            for machine_type in list(machine_types.keys()):
                if machine_type not in names:
                    to_remove.add(machine_type)

        for machine_type in to_remove:
            del machine_types[machine_type]
        return machine_types

    def _add_gpu(self, gpu_count):
        """Add a number of NVIDIA gpus to the current executor. This works
           by way of adding nvidia_gpu to the job default resources, and also
           changing the default machine type prefix to be n1, which is
           the currently only supported instance type for using GPUs for LHS.
        """
        if not gpu_count or gpu_count == 0:
            return

        logger.debug(
            "found resource request for {} GPUs. This will limit to n1 "
            "instance types.".format(gpu_count)
        )
        self.workflow.default_resources.parsed["nvidia_gpu"] = gpu_count
        self.workflow.default_resources.args.append("nvidia_gpu=%s" % gpu_count)
        self._machine_type_prefix = self._machine_type_prefix or ""
        if not self._machine_type_prefix.startswith("n1"):
            self._machine_type_prefix = "n1"

    def _generate_job_resources(self, job):
        """given a particular job, generate the resources that it needs,
           including default regions and the virtual machine configuration
        """
        # Right now, do a best effort mapping of resources to instance types
        cores = job.resources.get("_cores", 1)
        mem_mb = job.resources.get("mem_mb", 15360)
        disk_mb = job.resources.get("disk_mb", 128000)

        # Convert mb to gb, add buffer of 50
        disk_gb = math.ceil(disk_mb / 1024) + 10

        # Look for if the user wants an nvidia gpu
        gpu_count = job.resources.get("nvidia_gpu") or job.resources.get("gpu")

        # Update default resources using decided memory and disk
        self.workflow.default_resources = self.default_resources
        self.workflow.default_resources.args = [
            "mem_mb=%s" % mem_mb,
            "disk_mb=%s" % disk_mb,
        ]
        self.workflow.default_resources.parsed["mem_mb"] = mem_mb
        self.workflow.default_resources.parsed["disk_mb"] = disk_mb

        # If gpu wanted, limit to N1 general family, and update arguments
        if gpu_count:
            self._add_gpu(gpu_count)

        machine_types = self.get_available_machine_types()

        # Alert the user of machine_types available before filtering
        # https://cloud.google.com/compute/docs/machine-types
        logger.debug(
            "found {} machine types across regions {} before filtering "
            "to increase selection, define fewer regions".format(
                len(machine_types), self.regions
            )
        )

        # First pass - eliminate anything that too low in cpu/memory
        keepers = dict()

        # Also keep track of max cpus and memory, in case none available
        max_cpu = 1
        max_mem = 15360

        for name, machine_type in machine_types.items():
            max_cpu = max(max_cpu, machine_type["guestCpus"])
            max_mem = max(max_mem, machine_type["memoryMb"])
            if machine_type["guestCpus"] < cores or machine_type["memoryMb"] < mem_mb:
                continue
            keepers[name] = machine_type

        # If a prefix is set, filter down to it
        if self._machine_type_prefix:
            machine_types = keepers
            keepers = dict()
            for name, machine_type in machine_types.items():
                if name.startswith(self._machine_type_prefix):
                    keepers[name] = machine_type

        # If we don't have any contenders, workflow error
        if not keepers:
            if self._machine_type_prefix:
                raise WorkflowError(
                    "Machine prefix {prefix} is too strict, or the resources cannot "
                    " be satisfied, so there are no options "
                    "available.".format(prefix=self._machine_type_prefix)
                )
            else:
                raise WorkflowError(
                    "You requested {requestMemory} MB memory, {requestCpu} cores. "
                    "The maximum available are {availableMemory} MB memory and "
                    "{availableCpu} cores. These resources cannot be satisfied. "
                    "Please consider reducing the resource requirements of the "
                    "corresponding rule.".format(
                        requestMemory=mem_mb,
                        requestCpu=cores,
                        availableCpu=max_cpu,
                        availableMemory=max_mem,
                    )
                )

        # Now find (quasi) minimal to satisfy constraints
        machine_types = keepers

        # Select the first as the "smallest"
        smallest = list(machine_types.keys())[0]
        min_cores = machine_types[smallest]["guestCpus"]
        min_mem = machine_types[smallest]["memoryMb"]

        for name, machine_type in machine_types.items():
            if (
                machine_type["guestCpus"] < min_cores
                and machine_type["memoryMb"] < min_mem
            ):
                smallest = name
                min_cores = machine_type["guestCpus"]
                min_mem = machine_type["memoryMb"]

        selected = machine_types[smallest]
        logger.debug(
            "Selected machine type {}:{}".format(smallest, selected["description"])
        )
        virtual_machine = {
            "machineType": smallest,
            "labels": {"app": "snakemake"},
            "bootDiskSizeGb": disk_gb,
        }

        # If the user wants gpus, add accelerators here
        if gpu_count:
            accelerator = self._get_accelerator(gpu_count, zone=selected["zone"])
            virtual_machine["accelerators"] = [
                {"type": accelerator["name"], "count": gpu_count}
            ]

        resources = {"regions": self.regions, "virtualMachine": virtual_machine}
        return resources

    def _get_accelerator(self, gpu_count, zone):
        """Get an appropriate accelerator for a GPU given a zone selection.
           Currently Google offers NVIDIA Tesla T4 (likely the best),
           NVIDIA P100, and the same T4 for a graphical workstation. Since
           this isn't a graphical workstation use case, we choose the 
           accelerator that has >= to the maximumCardsPerInstace
        """
        if not gpu_count or gpu_count == 0:
            return

        accelerators = self._retry_request(
            self._compute_cli.acceleratorTypes().list(project=self.project, zone=zone)
        )

        # Filter down to those with greater than or equal to needed gpus
        keepers = {}
        for accelerator in accelerators.get("items", []):

            # We don't need a virtual workstation
            if accelerator["name"].endswith("vws"):
                continue

            if accelerator["maximumCardsPerInstance"] >= gpu_count:
                keepers[accelerator["name"]] = accelerator

        # Find smallest (in future the user might have preference for the type)
        smallest = list(keepers.keys())[0]
        max_gpu = keepers[smallest]["maximumCardsPerInstance"]

        # This should usually return P-100, which would be preference (cheapest)
        for name, accelerator in keepers.items():
            if accelerator["maximumCardsPerInstance"] < max_gpu:
                smallest = name
                max_gpu = accelerator["maximumCardsPerInstance"]

        return keepers[smallest]

    def _set_snakefile(self):
        """The snakefile must be a relative path, which cannot be reliably
           derived from the self.workflow.snakefile as we might have moved
           execution into a temporary directory, and the initial Snakefile
           was somewhere else on the system.
        """
        from snakemake import SNAKEFILE_CHOICES

        for snakefile in SNAKEFILE_CHOICES:
            if os.path.exists(os.path.join(self.workdir, snakefile)):
                self.snakefile = snakefile
                break

    def _set_workflow_sources(self, warning_size_gb=20):
        """We only add files from the working directory that are config related
           (e.g., the Snakefile or a config.yml equivalent), or checked into git. 
           Additionally, given that we encourage these packages to be small,
           we set a warning at 20GB.
        """

        def check_size(f):
            """A helper function to check the filesize, and return the file
               to the calling function
            """
            gb = bytesto(os.stat(f).st_size, "g")
            if gb > warning_size_gb:
                logger.warning(
                    "File {} (size {} GB) is greater than the {} GB suggested size "
                    "Consider uploading larger files to storage first.".format(
                        f, gb, warning_size_gb
                    )
                )
            return f

        self.workflow_sources = []

        for wfs in self.workflow.get_sources():
            if os.path.isdir(wfs):
                for (dirpath, dirnames, filenames) in os.walk(wfs):
                    self.workflow_sources.extend(
                        [check_size(os.path.join(dirpath, f)) for f in filenames]
                    )
            else:
                self.workflow_sources.append(check_size(os.path.abspath(wfs)))

    def _generate_build_source_package(self):
        """in order for the instance to access the working directory in storage,
           we need to upload it. This file is cleaned up at the end of the run.
           We do this, and then obtain from the instance and extract.
        """
        # Workflow sources for cloud executor must all be under same workdir root
        for filename in self.workflow_sources:
            if self.workdir not in filename:
                raise WorkflowError(
                    "All source files must be present in the working directory, "
                    "{workdir} to be uploaded to a build package that respects "
                    "relative paths, but {filename} was found outside of this "
                    "directory. Please set your working directory accordingly, "
                    "and the path of your Snakefile to be relative to it.".format(
                        workdir=self.workdir, filename=filename
                    )
                )

        # We will generate a tar.gz package, renamed by hash
        tmpname = next(tempfile._get_candidate_names())
        targz = os.path.join(tempfile.gettempdir(), "snakemake-%s.tar.gz" % tmpname)
        tar = tarfile.open(targz, "w:gz")

        # Add all workflow_sources files
        for filename in self.workflow_sources:
            arcname = filename.replace(self.workdir + os.path.sep, "")
            tar.add(filename, arcname=arcname)

        tar.close()

        # Rename based on hash, in case user wants to save cache
        sha256 = get_file_hash(targz)
        hash_tar = os.path.join(self.workdir, "snakeworkdir-%s.tar.gz" % sha256)

        # Only copy if we don't have it yet, clean up if we do
        if not os.path.exists(hash_tar):
            shutil.move(targz, hash_tar)
        else:
            os.remove(targz)

        # We will clean these all up at shutdown
        self._build_packages.add(hash_tar)

        return hash_tar

    def _generate_job_action(self, job):
        """generate a single action to execute the job.
        """
        # Derive the entrypoint command, the same content that might be written by self.get_jobscript(job)
        use_threads = "--force-use-threads" if not job.is_group() else ""

        exec_job = self.format_job(
            self.exec_job, job, _quote_all=True, use_threads=use_threads
        )

        # Now that we've parsed the job resource requirements, add to exec
        exec_job += self.get_default_resources_args()

        # The full command to download the archive, extract, and run
        # For snakemake bases, we must activate the conda environment, but
        # for custom images we must allow this to fail (hence || true)
        commands = [
            "/bin/bash",
            "-c",
            "mkdir -p /workdir && cd /workdir && wget -O /download.py https://gist.githubusercontent.com/vsoch/84886ef6469bedeeb9a79a4eb7aec0d1/raw/181499f8f17163dcb2f89822079938cbfbd258cc/download.py && chmod +x /download.py && source activate snakemake || true && pip install crc32c && python /download.py download %s %s /tmp/workdir.tar.gz && tar -xzvf /tmp/workdir.tar.gz && %s"
            % (self.bucket.name, self.pipeline_package, exec_job),
        ]

        # We are only generating one action, one job per run
        # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.pipelines/run#Action
        action = {
            "containerName": "snakejob-{}-{}".format(job.name, job.jobid),
            "imageUri": self.container_image,
            "commands": commands,
            "environment": self._generate_environment(),
            "labels": self._generate_pipeline_labels(job),
        }
        return action

    def _get_jobname(self, job):
        # Use a dummy job name (human readable and also namespaced)
        return "snakejob-%s-%s-%s" % (self.run_namespace, job.name, job.jobid)

    def _generate_pipeline_labels(self, job):
        """generate basic labels to identify the job, namespace, and that 
           snakemake is running the show!
        """
        jobname = self._get_jobname(job)
        labels = {"name": jobname, "app": "snakemake"}
        return labels

    def _generate_environment(self):
        """loop through envvars (keys to host environment) and add
           any that are requested for the container environment.
        """
        envvars = {}
        for key in self.envvars:
            try:
                envvars[key] = os.environ[key]
            except KeyError:
                continue

        # Warn the user that we cannot support secrets
        if envvars:
            logger.warning("This API does not support environment secrets.")
        return envvars

    def _generate_pipeline(self, job):
        """based on the job details, generate a google Pipeline object
           to pass to pipelines.run. This includes actions, resources,
           environment, and timeout.
        """
        # Generate targz with Snakefile and other context
        targz = self._generate_build_source_package()

        # Upload to temporary storage, only if doesn't exist
        self.pipeline_package = "source/cache/%s" % os.path.basename(targz)
        blob = self.bucket.blob(self.pipeline_package)
        logger.debug("build-package=%s" % self.pipeline_package)
        if not blob.exists():
            blob.upload_from_filename(targz, content_type="application/gzip")

        # Generate actions (one per job) and resources
        resources = self._generate_job_resources(job)
        action = self._generate_job_action(job)

        pipeline = {
            # Ordered list of actions to execute
            "actions": [action],
            # resources required for execution
            "resources": resources,
            # Technical question - difference between resource and action environment
            # For now we will set them to be the same.
            "environment": self._generate_environment(),
        }

        # "timeout": string in seconds (3.5s) is not included (defaults to 7 days)
        return pipeline

    def run(self, job, callback=None, submit_callback=None, error_callback=None):

        super()._run(job)

        # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.pipelines
        pipelines = self._api.projects().locations().pipelines()

        # pipelines.run
        # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.pipelines/run

        labels = self._generate_pipeline_labels(job)
        pipeline = self._generate_pipeline(job)

        # The body of the request is a Pipeline and labels
        body = {"pipeline": pipeline, "labels": labels}

        # capabilities - this won't currently work (Singularity in Docker)
        # We either need to add CAPS or run in privileged mode (ehh)
        if job.needs_singularity and self.workflow.use_singularity:
            logger.warning(
                "Singularity requires additional capabilities that "
                "aren't yet supported for standard Docker runs."
            )

        # location looks like: "projects/<project>/locations/<location>"
        operation = pipelines.run(parent=self.location, body=body)

        # 403 will result if no permission to use pipelines or project
        result = self._retry_request(operation)

        # The jobid is the last number of the full name
        jobid = result["name"].split("/")[-1]

        # Give some logging for how to get status
        logger.info(
            "Get status with:\n"
            "gcloud config set project {project}\n"
            "gcloud beta lifesciences operations describe {jobid}\n"
            "gcloud beta lifesciences operations list".format(
                project=self.project, jobid=jobid
            )
        )

        self.active_jobs.append(
            GoogleLifeSciencesJob(job, result["name"], jobid, callback, error_callback)
        )

    def _job_was_successful(self, status):
        """based on a status response (a [pipeline].projects.locations.operations.get
           debug print the list of events, return True if all return codes 0
           and False otherwise (indication of failure). In that a nonzero exit
           status is found, we also debug print it for the user.
        """
        success = True

        # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/Event
        for event in status["metadata"]["events"]:

            logger.debug(event["description"])

            # Does it always result in fail for other failure reasons?
            if "failed" in event:
                success = False
                action = event.get("failed")
                logger.error("{}: {}".format(action["code"], action["cause"]))

            elif "unexpectedExitStatus" in event:
                action = event.get("unexpectedExitStatus")

                if action["exitStatus"] != 0:
                    success = False

                    # Provide reason for the failure (desc includes exit code)
                    msg = "%s" % event["description"]
                    if "stderr" in action:
                        msg += ": %s" % action["stderr"]
                        logger.debug(msg)

        return success

    def _retry_request(self, request, timeout=2, attempts=3):
        """The Google Python API client frequently has BrokenPipe errors. This
           function takes a request, and executes it up to number of retry,
           each time with a 2* increase in timeout.

           Parameters
           ==========
           request: the Google Cloud request that needs to be executed
           timeout: time to sleep (in seconds) before trying again
           attempts: remaining attempts, throw error when hit 0
        """
        import googleapiclient

        try:
            return request.execute()
        except BrokenPipeError as ex:
            if attempts > 0:
                time.sleep(timeout)
                return self._retry_request(request, timeout * 2, attempts - 1)
            raise ex
        except googleapiclient.errors.HttpError as ex:
            log_verbose_traceback(ex)
            raise ex
        except (Exception, BaseException) as ex:
            log_verbose_traceback(ex)
            raise ex

    def _wait_for_jobs(self):
        """wait for jobs to complete. This means requesting their status,
           and then marking them as finished when a "done" parameter
           shows up. Even for finished jobs, the status should still return
        """
        import googleapiclient

        while True:
            # always use self.lock to avoid race conditions
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()

            # Loop through active jobs and act on status
            for j in active_jobs:

                # use self.status_rate_limiter to avoid too many API calls.
                with self.status_rate_limiter:

                    # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.operations/get
                    # Get status from projects.locations.operations/get
                    operations = self._api.projects().locations().operations()
                    request = operations.get(name=j.jobname)
                    logger.debug("Checking status for operation {}".format(j.jobid))

                    try:
                        status = self._retry_request(request)
                    except googleapiclient.errors.HttpError as ex:

                        # Operation name not found, even finished should be found
                        if ex.status == 404:
                            j.error_callback(j.job)
                            continue

                        # Unpredictable server (500) error
                        elif ex.status == 500:
                            logger.error(ex["content"].decode("utf-8"))
                            j.error_callback(j.job)

                    except WorkflowError as ex:
                        print_exception(ex, self.workflow.linemaps)
                        j.error_callback(j.job)
                        continue

                    # The operation is done
                    if status.get("done", False) == True:

                        # Derive success/failure from status codes (prints too)
                        if self._job_was_successful(status):
                            j.callback(j.job)
                        else:
                            self.print_job_error(j.job, jobid=j.jobid)
                            j.error_callback(j.job)

                    # The operation is still running
                    else:
                        still_running.append(j)

            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()
