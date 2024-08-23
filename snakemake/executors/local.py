__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


import os
import sys
import time
import shlex
import concurrent.futures
import subprocess
from functools import partial
from snakemake.common import async_run
from snakemake.executors import change_working_directory
from snakemake.settings.types import DeploymentMethod

from snakemake_interface_executor_plugins.executors.base import SubmittedJobInfo
from snakemake_interface_executor_plugins.executors.real import RealExecutor
from snakemake_interface_executor_plugins.dag import DAGExecutorInterface
from snakemake_interface_executor_plugins.workflow import WorkflowExecutorInterface
from snakemake_interface_executor_plugins.logging import LoggerExecutorInterface
from snakemake_interface_executor_plugins.jobs import (
    JobExecutorInterface,
    SingleJobExecutorInterface,
    GroupJobExecutorInterface,
)
from snakemake_interface_executor_plugins.settings import ExecMode, CommonSettings

from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.exceptions import print_exception, get_exception_origin
from snakemake.exceptions import format_error, RuleException, log_verbose_traceback
from snakemake.exceptions import (
    WorkflowError,
    SpawnedJobError,
    CacheMissException,
)


common_settings = CommonSettings(
    non_local_exec=False,
    implies_no_shared_fs=False,
    job_deploy_sources=False,
    pass_envvar_declarations_to_cmd=True,
    auto_deploy_default_storage_provider=False,
)


_ProcessPoolExceptions = (KeyboardInterrupt,)
try:
    from concurrent.futures.process import BrokenProcessPool

    _ProcessPoolExceptions = (KeyboardInterrupt, BrokenProcessPool)
except ImportError:
    pass


class Executor(RealExecutor):
    def __post_init__(self):
        self.use_threads = self.workflow.execution_settings.use_threads
        self.keepincomplete = self.workflow.execution_settings.keep_incomplete
        cores = self.workflow.resource_settings.cores

        # Zero thread jobs do not need a thread, but they occupy additional workers.
        # Hence we need to reserve additional workers for them.
        workers = cores + 5 if cores is not None else 5
        self.workers = workers
        self.pool = concurrent.futures.ThreadPoolExecutor(max_workers=self.workers)

    def get_exec_mode(self):
        return ExecMode.SUBPROCESS

    @property
    def job_specific_local_groupid(self):
        return False

    def get_job_exec_prefix(self, job: JobExecutorInterface):
        return f"cd {shlex.quote(self.workflow.workdir_init)}"

    def get_python_executable(self):
        return sys.executable

    def additional_general_args(self):
        return "--quiet progress rules"

    def get_job_args(self, job: JobExecutorInterface, **kwargs):
        return f"{super().get_job_args(job, **kwargs)}"

    def run_job(
        self,
        job: JobExecutorInterface,
    ):
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

        job_info = SubmittedJobInfo(job=job)

        future.add_done_callback(partial(self._callback, job_info))
        self.report_job_submission(job_info)

    def job_args_and_prepare(self, job: JobExecutorInterface):
        async_run(job.prepare())

        conda_env = (
            job.conda_env.address
            if DeploymentMethod.CONDA
            in self.workflow.deployment_settings.deployment_method
            and job.conda_env
            else None
        )
        container_img = (
            job.container_img_path
            if DeploymentMethod.APPTAINER
            in self.workflow.deployment_settings.deployment_method
            else None
        )
        env_modules = (
            job.env_modules
            if DeploymentMethod.ENV_MODULES
            in self.workflow.deployment_settings.deployment_method
            else None
        )

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
            self.workflow.output_settings.benchmark_extended,
            conda_env,
            container_img,
            self.workflow.deployment_settings.apptainer_args,
            env_modules,
            DeploymentMethod.APPTAINER
            in self.workflow.deployment_settings.deployment_method,
            self.workflow.linemaps,
            self.workflow.execution_settings.debug,
            self.workflow.execution_settings.cleanup_scripts,
            job.shadow_dir,
            job.jobid,
            (
                self.workflow.execution_settings.edit_notebook
                if self.dag.is_edit_notebook_job(job)
                else None
            ),
            self.workflow.conda_base_path,
            job.rule.basedir,
            self.workflow.sourcecache.cache_path,
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

        service_futures = [self.run_single_job(j) for j in job if j.is_service]
        normal_futures = [self.run_single_job(j) for j in job if not j.is_service]

        while normal_futures:
            for f in list(normal_futures):
                if f.done():
                    logger.debug("Job inside group is finished.")
                    ex = f.exception()
                    if ex is not None:
                        logger.debug(f"Job inside group failed with exception {ex}.")
                        # kill all shell commands of the other group jobs
                        # there can be only shell commands because the
                        # run directive is not allowed for pipe jobs
                        for j in job:
                            shell.kill(j.jobid)
                        raise ex
                    normal_futures.remove(f)
            time.sleep(1)

        if service_futures:
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

    def spawn_job(self, job: SingleJobExecutorInterface):
        cmd = self.format_job_exec(job)
        logger.debug(f"spawned job: {cmd}")

        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            raise SpawnedJobError()

    def cached_or_run(self, job: SingleJobExecutorInterface, run_func, *args):
        """
        Either retrieve result from cache, or run job with given function.
        """
        cache_mode = self.workflow.get_cache_mode(job.rule)
        try:
            if cache_mode:
                async_run(self.workflow.output_file_cache.fetch(job, cache_mode))
                return
        except CacheMissException:
            pass
        run_func(*args)
        if cache_mode:
            async_run(self.workflow.output_file_cache.store(job, cache_mode))

    def shutdown(self):
        self.pool.shutdown()

    def cancel(self):
        self.pool.shutdown()

    def _callback(self, job_info: SubmittedJobInfo, future):
        try:
            ex = future.exception()
            if ex is not None:
                raise ex
            self.report_job_success(job_info)
        except _ProcessPoolExceptions:
            self.handle_job_error(job_info.job)
            # no error callback, just silently ignore the interrupt as the main scheduler is also killed
        except SpawnedJobError:
            # don't print error message, this is done by the spawned subprocess
            self.report_job_error(job_info)
        except Exception as ex:
            if self.workflow.output_settings.verbose or (
                not job_info.job.is_group() and not job_info.job.is_shell
            ):
                print_exception(ex, self.workflow.linemaps)
            self.report_job_error(job_info)

    @property
    def cores(self):
        return self.workflow.resource_settings.cores


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
    benchmark_extended,
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
    sourcecache_path,
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
                            sourcecache_path,
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
                                sourcecache_path,
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
                    sourcecache_path,
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
            # Add job info to (all repeats of) benchmark file
            for bench_record in bench_records:
                bench_record.jobid = jobid
                bench_record.rule_name = job_rule.name
                bench_record.wildcards = wildcards
                bench_record.params = params
                bench_record.resources = resources
                bench_record.input = input
                bench_record.threads = threads
            write_benchmark_records(bench_records, benchmark, benchmark_extended)
        except Exception as ex:
            raise WorkflowError(ex)
