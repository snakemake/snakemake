from itertools import chain
from snakemake.benchmark import BenchmarkRecord
from snakemake.shell import shell
from typing import Any
from typing import List
from typing import Union
from typing import Dict
from typing import Optional

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
from snakemake.io import Namedlist
from attr import dataclass
import os
import sys
import time
import shlex
import concurrent.futures
import subprocess
from functools import partial
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
from snakemake_interface_software_deployment_plugins import EnvBase as SoftwareEnvBase

from snakemake.rules import Rule
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
        return "--quiet progress rules host"

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

    def job_args_and_prepare(self, job: JobExecutorInterface) -> "RunArgs":
        self.workflow.async_run(job.prepare())

        benchmark = None
        benchmark_repeats = job.benchmark_repeats or 1
        if job.benchmark is not None:
            benchmark = str(job.benchmark)

        return RunArgs(
            job_rule=job.rule,
            input=job.input._plainstrings(),
            output=job.output._plainstrings(),
            params=job.params,
            wildcards=job.wildcards,
            threads=job.threads,
            resources=job.resources,
            log=job.log._plainstrings(),
            benchmark=benchmark,
            benchmark_repeats=benchmark_repeats,
            benchmark_extended=self.workflow.output_settings.benchmark_extended,
            linemaps=self.workflow.linemaps,
            debug=self.workflow.execution_settings.debug,
            cleanup_scripts=self.workflow.execution_settings.cleanup_scripts,
            shadow_dir=job.shadow_dir,
            jobid=job.jobid,
            edit_notebook=(
                self.workflow.execution_settings.edit_notebook
                if self.dag.is_edit_notebook_job(job)
                else None
            ),
            basedir=job.rule.basedir,
            cache_path=self.workflow.sourcecache.cache_path,
            runtime_cache_path=self.workflow.sourcecache.runtime_cache_path,
            runtime_paths=self.workflow.runtime_paths,
            wrapper_prefix=self.workflow.workflow_settings.wrapper_prefix,
            not_block_search_path_envvars=(
                self.workflow.deployment_settings.not_block_search_path_envvars
            ),
            software_env=job.software_env,
        )

    def run_single_job(self, job: SingleJobExecutorInterface):
        if (
            self.use_threads
            or (not job.is_shadow and not job.is_run)
            or job.is_template_engine
        ):
            future = self.pool.submit(
                self.cached_or_run, job, run_wrapper, self.job_args_and_prepare(job)
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
        self.workflow.async_run(self.acached_or_run(job, run_func, *args))

    async def acached_or_run(self, job: SingleJobExecutorInterface, run_func, *args):
        """
        Either retrieve result from cache, or run job with given function.
        """
        cache_mode = self.workflow.get_cache_mode(job.rule)
        try:
            if cache_mode:
                await self.workflow.output_file_cache.fetch(job, cache_mode)
                return
        except CacheMissException:
            pass
        run_func(*args)
        if cache_mode:
            await self.workflow.output_file_cache.store(job, cache_mode)

    def shutdown(self):
        self.pool.shutdown()

    def cancel(self):
        self.pool.shutdown()

    def _callback(self, job_info: SubmittedJobInfo, future):
        try:
            ex = future.exception()
            if ex is not None:
                print_exception(ex, self.workflow.linemaps)
                self.report_job_error(job_info)
            else:
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


@dataclass
class RunArgs:
    job_rule: Rule
    input: Namedlist
    output: Namedlist
    params: Namedlist
    wildcards: Namedlist
    threads: int
    resources: Namedlist
    log: Namedlist
    benchmark: Optional[str]
    benchmark_repeats: int
    benchmark_extended: bool
    linemaps: Dict[str, Dict[int, int]]
    debug: bool
    cleanup_scripts: bool
    shadow_dir: Path
    jobid: Union[int, str]
    edit_notebook: bool
    basedir: Path
    cache_path: Path
    runtime_cache_path: Path
    runtime_paths: List[Path]
    not_block_search_path_envvars: bool
    wrapper_prefix: Optional[str]
    software_env: Optional[SoftwareEnvBase]
    bench_iteration: int = 0
    bench_record: Optional[BenchmarkRecord] = None

    def rulefunc_args(self) -> Dict[str, Any]:
        return {
            item: getattr(self, item) for item in self.rulefunc_args_attributes()
        } | {"rule": self.job_rule.name, "run_args": self}

    @classmethod
    def rulefunc_args_attributes(cls) -> List[str]:
        return [
            "input",
            "output",
            "params",
            "wildcards",
            "threads",
            "resources",
            "log",
            "jobid",
            "bench_iteration",
        ]

    @classmethod
    def rulefunc_args_signature(cls) -> str:
        return ", ".join(chain(cls.rulefunc_args_attributes(), ["rule", "run_args"]))


    @property
    def is_shell(self) -> bool:
        return self.job_rule.shellcmd is not None


def run_wrapper(run_args: RunArgs):
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
    run = run_args.job_rule.run_func

    if os.name == "posix" and run_args.debug:
        sys.stdin = open("/dev/stdin")

    if run_args.benchmark is not None:
        from snakemake.benchmark import (
            benchmarked,
            write_benchmark_records,
        )

    # Change workdir if shadow defined and not using container.
    # Otherwise, we do the change from inside the container.
    immediate_shadow_dir = None
    if run_args.software_env is None or run_args.software_env.spec.kind != "container":
        immediate_shadow_dir = run_args.shadow_dir
        run_args.shadow_dir = None

    try:
        with change_working_directory(immediate_shadow_dir):

            if run_args.benchmark:
                bench_records = []
                for bench_iteration in range(run_args.benchmark_repeats):
                    run_args.bench_iteration = bench_iteration
                    # Determine whether to benchmark this process or do not
                    # benchmarking at all.  We benchmark this process unless the
                    # execution is done through the ``shell:``, ``script:``, or
                    # ``wrapper:`` stanza.
                    is_sub = (
                        run_args.job_rule.shellcmd
                        or run_args.job_rule.script
                        or run_args.job_rule.wrapper
                        or run_args.job_rule.cwl
                    )
                    if is_sub:
                        # The benchmarking through ``benchmarked()`` is started
                        # in the execution of the shell fragment, script, wrapper
                        # etc, as the child PID is available there.
                        run_args.bench_record = BenchmarkRecord()
                        run(**run_args.rulefunc_args())
                    else:
                        # The benchmarking is started here as we have a run section
                        # and the generated Python function is executed in this
                        # process' thread.
                        with benchmarked() as bench_record:
                            run_args.bench_record = bench_record
                            run(**run_args.rulefunc_args())
                    # Store benchmark record for this iteration
                    bench_records.append(bench_record)
            else:
                run(**run_args.rulefunc_args())
    except (KeyboardInterrupt, SystemExit) as e:
        # Re-raise the keyboard interrupt in order to record an error in the
        # scheduler but ignore it
        raise e
    except BaseException as ex:
        # this ensures that exception can be re-raised in the parent thread
        origin = get_exception_origin(ex, run_args.linemaps)
        if origin is not None:
            log_verbose_traceback(ex)
            lineno, file = origin
            raise RuleException(
                format_error(
                    ex,
                    lineno,
                    linemaps=run_args.linemaps,
                    snakefile=file,
                    show_traceback=True,
                )
            ) from ex
        else:
            # some internal bug, just reraise
            raise ex

    if run_args.benchmark is not None:
        try:
            # Add job info to (all repeats of) benchmark file
            for bench_record in bench_records:
                bench_record.jobid = run_args.jobid
                bench_record.rule_name = run_args.job_rule.name
                bench_record.wildcards = run_args.wildcards
                bench_record.params = run_args.params
                bench_record.resources = run_args.resources
                bench_record.input = run_args.input
                bench_record.threads = run_args.threads
            write_benchmark_records(
                bench_records, run_args.benchmark, run_args.benchmark_extended
            )
        except Exception as ex:
            raise WorkflowError(ex) from ex
