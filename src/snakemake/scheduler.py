__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from bisect import bisect
from collections import defaultdict, deque
import math
import os, signal, sys
import threading

from itertools import chain, accumulate, repeat
from contextlib import ContextDecorator
import time

from snakemake_interface_executor_plugins.scheduler import JobSchedulerExecutorInterface
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
from snakemake_interface_executor_plugins.registry import Plugin as ExecutorPlugin
from snakemake_interface_executor_plugins.settings import ExecMode
from snakemake_interface_logger_plugins.common import LogEvent
from snakemake.common import async_run

from snakemake.exceptions import RuleException, WorkflowError, print_exception
from snakemake.settings.enums import Quietness
from snakemake.logging import logger
from snakemake.jobs import GroupJob

from snakemake.settings.types import MaxJobsPerTimespan

registry = ExecutorPluginRegistry()


def cumsum(iterable, zero=[0]):
    return list(chain(zero, accumulate(iterable)))


_ERROR_MSG_FINAL = (
    "Exiting because a job execution failed. Look below for error messages"
)

_ERROR_MSG_ISSUE_823 = (
    "BUG: Out of jobs ready to be started, but not all files built yet."
    " Please check https://github.com/snakemake/snakemake/issues/823 for more information."
)


class DummyRateLimiter(ContextDecorator):
    async def __aenter__(self):
        return self

    async def __aexit__(self, *args):
        return False


class JobScheduler(JobSchedulerExecutorInterface):
    def __init__(self, workflow, executor_plugin: ExecutorPlugin):
        """Create a new instance of KnapsackJobScheduler."""
        self.workflow = workflow

        self.dryrun = self.workflow.dryrun
        self.touch = self.workflow.touch
        self.quiet = self.workflow.output_settings.quiet
        self.keepgoing = self.workflow.execution_settings.keep_going
        self.running = set()
        self.failed = set()
        self.finished_jobs = 0
        self.greediness = self.workflow.scheduling_settings.greediness
        self.subsample = self.workflow.scheduling_settings.subsample
        self._tofinish = []
        self._toerror = []
        self.handle_job_success = True
        self.update_resources = True
        self.print_progress = (
            not self.quiet or Quietness.PROGRESS not in self.quiet
        ) and not self.dryrun
        self.update_checkpoint_dependencies = not self.dryrun
        self.job_rate_limiter = (
            JobRateLimiter(self.workflow.scheduling_settings.max_jobs_per_timespan)
            if not self.dryrun
            and self.workflow.scheduling_settings.max_jobs_per_timespan
            else None
        )

        nodes_unset = workflow.global_resources["_nodes"] is None

        self.global_resources = {
            name: (sys.maxsize if res is None else res)
            for name, res in workflow.global_resources.items()
        }

        if not nodes_unset:
            # Do not restrict cores locally if nodes are used (i.e. in case of cluster/cloud submission).
            self.global_resources["_cores"] = sys.maxsize
        # register job count resource (always initially unrestricted)
        self.global_resources["_job_count"] = sys.maxsize

        self.resources = dict(self.global_resources)

        self._open_jobs = threading.Semaphore(0)
        self._lock = threading.Lock()

        self._errors = False
        self._executor_error = None
        self._finished = False
        self._job_queue = None
        self._last_job_selection_empty = False
        self._last_update_queue_input_jobs = 0
        self.submit_callback = self._noop
        self.finish_callback = self._proceed
        self._run_performed = None

        if workflow.remote_execution_settings.immediate_submit:
            self.submit_callback = self._proceed
            self.finish_callback = self._noop

        self._local_executor = None

        if self.workflow.local_exec:
            self._executor = executor_plugin.executor(
                self.workflow,
                logger,
            )
        else:
            self._executor = executor_plugin.executor(
                self.workflow,
                logger,
            )
            self._local_executor = (
                ExecutorPluginRegistry()
                .get_plugin("local")
                .executor(
                    self.workflow,
                    logger,
                )
            )

        # Choose job selector (greedy or ILP)
        self._job_selector = self.job_selector_greedy
        if self.workflow.scheduling_settings.scheduler == "ilp":
            import pulp

            if pulp.apis.LpSolverDefault is None:
                logger.warning(
                    "Falling back to greedy scheduler because no default "
                    "solver is found for pulp (you have to install either "
                    "coincbc or glpk)."
                )
            else:
                self._job_selector = self.job_selector_ilp

        self._user_kill = None
        try:
            signal.signal(signal.SIGTERM, self.exit_gracefully)
        except ValueError:
            # If this fails, it is due to scheduler not being invoked in the main thread.
            # This can only happen with --gui, in which case it is fine for now.
            pass
        self._open_jobs.release()

    def executor_error_callback(self, exception):
        with self._lock:
            self._executor_error = exception
            # next scheduling round to catch and raise error
            self._open_jobs.release()

    @property
    def stats(self):
        return self._stats

    @property
    def open_jobs(self):
        """Return open jobs."""
        jobs = list(self.workflow.dag.ready_jobs)
        return jobs

    @property
    def remaining_jobs(self):
        """Return jobs to be scheduled including not yet ready ones."""
        return [
            job
            for job in self.workflow.dag.needrun_jobs()
            if job not in self.running
            and not self.workflow.dag.finished(job)
            and job not in self.failed
        ]

    def schedule(self):
        """Schedule jobs that are ready, maximizing cpu usage."""
        try:
            while True:
                if self.workflow.dag.queue_input_jobs:
                    self.update_queue_input_jobs()
                # work around so that the wait does not prevent keyboard interrupts
                # while not self._open_jobs.acquire(False):
                #    time.sleep(1)
                self._open_jobs.acquire()

                # obtain needrun and running jobs in a thread-safe way
                with self._lock:
                    self._finish_jobs()
                    self._error_jobs()
                    needrun = set(self.open_jobs)
                    running = list(self.running)
                    errors = self._errors
                    executor_error = self._executor_error
                    user_kill = self._user_kill

                # handle errors
                if user_kill or (not self.keepgoing and errors) or executor_error:
                    if user_kill == "graceful":
                        logger.info(
                            "Will exit after finishing currently running jobs (scheduler)."
                        )

                    if executor_error:
                        print_exception(executor_error, self.workflow.linemaps)

                    if executor_error or not running:
                        logger.info("Shutting down, this might take some time.")
                        self._executor.shutdown()
                        if not user_kill:
                            logger.error(_ERROR_MSG_FINAL)
                            for job in self.failed:
                                job.log_error()
                        return False
                    continue

                # all runnable jobs have finished, normal shutdown
                if (
                    not needrun
                    and (
                        not running
                        or self.workflow.remote_execution_settings.immediate_submit
                    )
                    and not self.workflow.dag.has_unfinished_queue_input_jobs()
                ):
                    self._executor.shutdown()
                    if errors:
                        logger.error(_ERROR_MSG_FINAL)
                        for job in self.failed:
                            job.log_error()
                    # we still have unfinished jobs. this is not good. direct
                    # user to github issue
                    if self.remaining_jobs and not self.keepgoing:
                        logger.error(_ERROR_MSG_ISSUE_823)
                        logger.error(
                            "Remaining jobs:\n"
                            + "\n".join(
                                " - " + str(job) + ": " + ", ".join(job.output)
                                for job in self.remaining_jobs
                            )
                        )
                        return False
                    return not errors

                # continue if no new job needs to be executed
                if not needrun:
                    if self.workflow.dag.has_unfinished_queue_input_jobs():
                        logger.info("Waiting for queue input...")
                        # schedule a reevaluation in 10 seconds
                        self._schedule_reevalutation(
                            self.workflow.execution_settings.queue_input_wait_time
                        )
                    continue

                # select jobs by solving knapsack problem (omit with dryrun)
                if self.dryrun:
                    run = needrun
                else:
                    # Reset params and resources because they might still contain TBDs
                    # or old values from before files have been regenerated.
                    # Now, they can be recalculated as all input is present and up to date.
                    for job in needrun:
                        job.reset_params_and_resources()

                    logger.debug(f"Resources before job selection: {self.resources}")

                    # Subsample jobs to be run (to speedup solver)
                    n_total_needrun = len(needrun)
                    if self.subsample and n_total_needrun > self.subsample:
                        import random

                        needrun = set(random.sample(tuple(needrun), k=self.subsample))
                        logger.debug(
                            f"Ready subsampled jobs: {len(needrun)} (out of {n_total_needrun})"
                        )
                    else:
                        logger.debug(f"Ready jobs: {n_total_needrun}")

                    if not self._last_job_selection_empty:
                        logger.info("Select jobs to execute...")
                    run = self.job_selector(needrun)
                    self._last_job_selection_empty = not run

                    logger.debug(f"Selected jobs: {len(run)}")
                    logger.debug(f"Resources after job selection: {self.resources}")

                # update running jobs
                with self._lock:
                    self.running.update(run)
                    # remove from ready_jobs
                    self.workflow.dag.register_running(run)

                if run:
                    if not self.dryrun:
                        logger.info(
                            f"Execute {len(run)} jobs...",
                            extra=dict(
                                event=LogEvent.JOB_STARTED, jobs=[j.jobid for j in run]
                            ),
                        )

                    # actually run jobs
                    local_runjobs = [job for job in run if job.is_local]
                    runjobs = [job for job in run if not job.is_local]
                    if local_runjobs:
                        if (
                            not self.workflow.remote_exec
                            and not self.workflow.local_exec
                        ):
                            # Workflow uses a remote plugin and this scheduling run
                            # is on the main process. Hence, we have to download
                            # non-shared remote files for the local jobs.
                            async_run(
                                self.workflow.dag.retrieve_storage_inputs(
                                    jobs=local_runjobs, also_missing_internal=True
                                )
                            )
                        self.run(
                            local_runjobs,
                            executor=self._local_executor or self._executor,
                        )
                    if runjobs:
                        self.run(runjobs)
                if not self.dryrun:

                    if self._run_performed is None or self._run_performed:
                        if self.running:
                            logger.debug("Waiting for running jobs to complete.")
                        else:
                            logger.debug("Waiting for more resources.")
                    if self.job_rate_limiter is not None:
                        # need to reevaluate because after the timespan we can
                        # schedule more jobs again
                        self._schedule_reevalutation(self.job_rate_limiter.timespan)
        except (KeyboardInterrupt, SystemExit):
            logger.info(
                "Terminating processes on user request, this might take some time."
            )
            self._executor.cancel()
            return False

    def _schedule_reevalutation(self, delay: int) -> None:
        threading.Timer(
            delay,
            lambda: self._open_jobs.release(),
        ).start()

    def _finish_jobs(self):
        # must be called from within lock
        # clear the global tofinish such that parallel calls do not interfere
        async def postprocess():
            for job in self._tofinish:
                # IMPORTANT: inside of this loop, there may be no calls that have
                # a complexity of at least the number of jobs.
                # Otherwise the function would be quadratic in the number of jobs.
                if not self.workflow.dryrun:
                    try:
                        if self.workflow.exec_mode == ExecMode.DEFAULT:
                            await job.postprocess(
                                store_in_storage=not self.touch,
                                handle_log=True,
                                handle_touch=not self.touch,
                                ignore_missing_output=self.touch,
                            )
                        elif self.workflow.exec_mode == ExecMode.SUBPROCESS:
                            await job.postprocess(
                                store_in_storage=False,
                                handle_log=True,
                                handle_touch=True,
                            )
                        else:
                            await job.postprocess(
                                # storage upload will be done after all jobs of
                                # this remote job (e.g. in case of group) are finished
                                # DAG.store_storage_outputs()
                                store_in_storage=False,
                                handle_log=True,
                                handle_touch=True,
                            )
                    except (RuleException, WorkflowError) as e:
                        # if an error occurs while processing job output,
                        # we do the same as in case of errors during execution
                        print_exception(e, self.workflow.linemaps)
                        await job.postprocess(error=True)
                        self._handle_error(job, postprocess_job=False)
                        continue

                if self.handle_job_success:
                    self.get_executor(job).handle_job_success(job)

                if self.update_resources:
                    # normal jobs have len=1, group jobs have len>1
                    self.finished_jobs += len(job)
                    self.running.remove(job)
                    self._free_resources(job)

                if self.print_progress:
                    if job.is_group():
                        for j in job:
                            logger.info(
                                f"Finished jobid: {j.jobid} (Rule: {j.rule.name})",
                                extra=dict(event=LogEvent.JOB_FINISHED, job_id=j.jobid),
                            )
                    else:
                        logger.info(
                            f"Finished jobid: {job.jobid} (Rule: {job.rule.name})",
                            extra=dict(event=LogEvent.JOB_FINISHED, job_id=job.jobid),
                        )
                    self.progress()

                await self.workflow.dag.finish(
                    job,
                    update_checkpoint_dependencies=self.update_checkpoint_dependencies,
                )

        async_run(postprocess())
        self._tofinish.clear()

    def update_queue_input_jobs(self):
        currtime = time.time()
        if currtime - self._last_update_queue_input_jobs >= 10:
            self._last_update_queue_input_jobs = currtime
            async_run(self.workflow.dag.update_queue_input_jobs())

    def _error_jobs(self):
        # must be called from within lock
        for job in self._toerror:
            self._handle_error(job)
        self._toerror.clear()

    def run(self, jobs, executor=None):
        self._run_performed = True
        if executor is None:
            executor = self._executor
        executor.run_jobs(jobs)

    def get_executor(self, job):
        if job.is_local and self._local_executor is not None:
            return self._local_executor
        return self._executor

    def _noop(self, job):
        pass

    def _free_resources(self, job):
        for name, value in job.scheduler_resources.items():
            if name in self.resources and name != "_job_count":
                value = self.calc_resource(name, value)
                self.resources[name] += value

    def _proceed(self, job):
        """Do stuff after job is finished."""
        with self._lock:
            self._tofinish.append(job)

            if self.dryrun:
                if len(self.running) - len(self._tofinish) - len(self._toerror) <= 0:
                    # During dryrun, only release when all running jobs are done.
                    # This saves a lot of time, as self.open_jobs has to be
                    # evaluated less frequently.
                    self._open_jobs.release()
            else:
                # go on scheduling if there is any free core
                self._open_jobs.release()

    def error_callback(self, job):
        with self._lock:
            self._toerror.append(job)
            self._open_jobs.release()

    def _handle_error(self, job, postprocess_job: bool = True):
        """Clear jobs and stop the workflow.

        If Snakemake is configured to restart jobs then the job might have
        "restart_times" left and we just decrement and let the scheduler
        try to run the job again.
        """
        # must be called from within lock
        if postprocess_job and not self.workflow.dryrun:
            async_run(
                job.postprocess(
                    error=True,
                )
            )
        self.get_executor(job).handle_job_error(job)
        self.running.remove(job)
        self._free_resources(job)
        # attempt starts counting from 1, but the first attempt is not
        # a restart, hence we subtract 1.
        if job.restart_times > job.attempt - 1:
            logger.info(f"Trying to restart job {self.workflow.dag.jobid(job)}.")
            job.attempt += 1
            # add job to those being ready again
            self.workflow.dag._ready_jobs.add(job)
        else:
            self._errors = True
            self.failed.add(job)

    def exit_gracefully(self, *args):
        with self._lock:
            self._user_kill = "graceful"
        self._open_jobs.release()

    def job_selector(self, jobs):
        # get number of free jobs to submit
        if self.job_rate_limiter is None:
            # ensure that the job count is not restricted
            assert (
                self.resources["_job_count"] == sys.maxsize
            ), f"Job count is {self.resources['_job_count']}, but should be {sys.maxsize}"
            return self._job_selector(jobs)
        n_free_jobs = self.job_rate_limiter.get_free_jobs()
        if n_free_jobs == 0:
            logger.info("Job rate limit reached, waiting for free slots.")
            return set()
        else:
            self.resources["_job_count"] = n_free_jobs
            selected = self._job_selector(jobs)
            # update job rate limiter
            self.job_rate_limiter.register_jobs(len(selected))
            return selected

    def job_selector_ilp(self, jobs):
        """
        Job scheduling by optimization of resource usage by solving ILP using pulp
        """
        import pulp
        from pulp import lpSum

        if len(jobs) == 1:
            return self.job_selector_greedy(jobs)

        with self._lock:
            if not self.resources["_cores"]:
                return set()

            scheduled_jobs = {
                job: pulp.LpVariable(
                    f"job_{idx}", lowBound=0, upBound=1, cat=pulp.LpInteger
                )
                for idx, job in enumerate(jobs)
            }

            temp_files = {
                temp_file
                for job in jobs
                for temp_file in self.workflow.dag.temp_input(job)
            }

            async def get_temp_sizes_gb():
                return {f: (await f.size()) / 1e9 for f in temp_files}

            temp_sizes_gb = (
                defaultdict(int) if self.touch else async_run(get_temp_sizes_gb())
            )

            temp_job_improvement = {
                temp_file: pulp.LpVariable(
                    f"temp_file_{idx}", lowBound=0, upBound=1, cat="Continuous"
                )
                for idx, temp_file in enumerate(temp_files)
            }

            temp_file_deletable = {
                temp_file: pulp.LpVariable(
                    f"deletable_{idx}",
                    lowBound=0,
                    upBound=1,
                    cat=pulp.LpInteger,
                )
                for idx, temp_file in enumerate(temp_files)
            }
            prob = pulp.LpProblem("JobScheduler", pulp.LpMaximize)

            total_temp_size = max(
                sum([temp_sizes_gb[temp_file] for temp_file in temp_files]), 1
            )
            total_core_requirement = sum(
                [max(job.scheduler_resources.get("_cores", 1), 1) for job in jobs]
            )
            # Objective function
            # Job priority > Core load
            # Core load > temp file removal
            # Instant removal > temp size
            prob += (
                2
                * total_core_requirement
                * 2
                * total_temp_size
                * lpSum([job.priority * scheduled_jobs[job] for job in jobs])
                + 2
                * total_temp_size
                * lpSum(
                    [
                        max(job.scheduler_resources.get("_cores", 1), 1)
                        * scheduled_jobs[job]
                        for job in jobs
                    ]
                )
                + total_temp_size
                * lpSum(
                    [
                        temp_file_deletable[temp_file] * temp_sizes_gb[temp_file]
                        for temp_file in temp_files
                    ]
                )
                + lpSum(
                    [
                        temp_job_improvement[temp_file] * temp_sizes_gb[temp_file]
                        for temp_file in temp_files
                    ]
                )
            )

            # Constraints:
            for name in self.global_resources:
                prob += (
                    lpSum(
                        [
                            scheduled_jobs[job] * job.scheduler_resources.get(name, 0)
                            for job in jobs
                        ]
                    )
                    <= self.resources[name]
                )

            # Choose jobs that lead to "fastest" (minimum steps) removal of existing temp file
            remaining_jobs = self.remaining_jobs
            for temp_file in temp_files:
                prob += temp_job_improvement[temp_file] <= lpSum(
                    [
                        scheduled_jobs[job] * self.required_by_job(temp_file, job)
                        for job in jobs
                    ]
                ) / lpSum(
                    [self.required_by_job(temp_file, job) for job in remaining_jobs]
                )

                prob += (
                    temp_file_deletable[temp_file] <= temp_job_improvement[temp_file]
                )

        status = self._solve_ilp(prob, time_limit=10)
        if pulp.LpStatus[status] != "Optimal":
            if pulp.LpStatus[status] == "Not Solved":
                logger.warning(
                    "Failed to solve scheduling problem with ILP solver in time (10s)."
                )
            elif pulp.LpStatus[status] == "Infeasible":
                logger.warning("Failed to solve scheduling problem with ILP solver.")

            logger.warning("Falling back to greedy solver.")
            return self.job_selector_greedy(jobs)

        selected_jobs = set(
            job
            for job, variable in scheduled_jobs.items()
            if math.isclose(variable.value(), 1.0)
        )

        if not selected_jobs:
            # No selected jobs. This could be due to insufficient resources or a failure in the ILP solver
            # Hence, we silently fall back to the greedy solver to make sure that we don't miss anything.
            return self.job_selector_greedy(jobs)

        self.update_available_resources(selected_jobs)
        return selected_jobs

    def _solve_ilp(self, prob, threads=2, time_limit=10):
        import pulp

        old_path = os.environ["PATH"]
        if self.workflow.scheduling_settings.solver_path is None:
            # Temporarily prepend the given snakemake env to the path, such that the solver can be found in any case.
            # This is needed for cluster envs, where the cluster job might have a different environment but
            # still needs access to the solver binary.
            os.environ["PATH"] = "{}:{}".format(
                self.workflow.scheduling_settings.solver_path,
                os.environ["PATH"],
            )
        try:
            solver = (
                pulp.getSolver(self.workflow.scheduling_settings.ilp_solver)
                if self.workflow.scheduling_settings.ilp_solver
                else pulp.apis.LpSolverDefault
            )
        finally:
            os.environ["PATH"] = old_path
        solver.optionsDict["threads"] = threads
        solver.timeLimit = time_limit
        solver.msg = self.workflow.output_settings.verbose
        return prob.solve(solver)

    def required_by_job(self, temp_file, job):
        return 1 if temp_file in self.workflow.dag.temp_input(job) else 0

    def job_selector_greedy(self, jobs):
        """
        Using the greedy heuristic from
        "A Greedy Algorithm for the General Multidimensional Knapsack
        Problem", Akcay, Li, Xu, Annals of Operations Research, 2012

        Args:
            jobs (list):    list of jobs
        """
        logger.debug("Selecting jobs to run using greedy solver.")

        with self._lock:
            if not self.resources["_cores"]:
                return set()
            # each job is an item with one copy (0-1 MDKP)
            n = len(jobs)
            x = [0] * n  # selected jobs
            E = set(range(n))  # jobs still free to select
            u = [1] * n
            a = list(map(self.job_weight, jobs))  # resource usage of jobs

            async def rewards():
                return [await self.job_reward(job) for job in jobs]

            c = async_run(rewards())  # job rewards

            def calc_reward():
                return [c_j * y_j for c_j, y_j in zip(c, y)]

            b = [
                self.resources[name] for name in self.global_resources
            ]  # resource capacities

            while True:
                # Step 2: compute effective capacities
                y = [
                    (
                        min(
                            (min(u[j], b_i // a_j_i) if a_j_i > 0 else u[j])
                            for b_i, a_j_i in zip(b, a[j])
                            if a_j_i
                        )
                        if j in E
                        else 0
                    )
                    for j in range(n)
                ]
                if not any(y):
                    break
                y = [
                    (max(1, int(self.greediness * y_j)) if y_j > 0 else 0) for y_j in y
                ]

                # Step 3: compute rewards on cumulative sums
                reward = calc_reward()
                j_sel = max(E, key=reward.__getitem__)  # argmax

                # Step 4: batch increment
                y_sel = y[j_sel]

                # Step 5: update information
                x[j_sel] += y_sel
                b = [b_i - (a_j_i * y_sel) for b_i, a_j_i in zip(b, a[j_sel])]
                u[j_sel] -= y_sel
                if not u[j_sel] or self.greediness == 1:
                    E.remove(j_sel)
                if not E:
                    break

            solution = set(job for job, sel in zip(jobs, x) if sel)

            self.update_available_resources(solution)
            return solution

    def update_available_resources(self, selected_jobs):
        for name in self.global_resources:
            # _job_count is updated per JobRateLimiter before scheduling
            if name != "_job_count":
                self.resources[name] -= sum(
                    [job.scheduler_resources.get(name, 0) for job in selected_jobs]
                )

    def calc_resource(self, name, value):
        gres = self.global_resources[name]
        if value > gres:
            if name == "_cores":
                name = "threads"
            raise WorkflowError(
                "Job needs {name}={res} but only {name}={gres} "
                "are available. This is likely because two "
                "jobs are connected via a pipe or a service output and have to run "
                "simultaneously. Consider providing more "
                "resources (e.g. via --cores).".format(name=name, res=value, gres=gres)
            )
        return value

    def rule_weight(self, rule):
        res = rule.resources
        return [
            self.calc_resource(name, res.get(name, 0)) for name in self.global_resources
        ]

    def job_weight(self, job):
        res = job.scheduler_resources
        return [
            self.calc_resource(name, res.get(name, 0)) for name in self.global_resources
        ]

    async def job_reward(self, job):
        if (
            self.touch
            or self.dryrun
            or self.workflow.remote_execution_settings.immediate_submit
        ):
            temp_size = 0
            input_size = 0
        else:
            try:
                temp_size = await self.workflow.dag.temp_size(job)
                input_size = await job.inputsize()
            except FileNotFoundError:
                # If the file is not yet present, this shall not affect the
                # job selection.
                temp_size = 0
                input_size = 0

        # Usually, this should guide the scheduler to first schedule all jobs
        # that remove the largest temp file, then the second largest and so on.
        # Since the weight is summed up, it can in theory be that it sometimes
        # prefers a set of many jobs that all depend on smaller temp files though.
        # A real solution to the problem is therefore to use dummy jobs that
        # ensure selection of groups of jobs that together delete the same temp
        # file.

        return (job.priority, temp_size, input_size)

    def progress(self):
        """Display the progress."""
        logger.info(
            None,
            extra=dict(
                event=LogEvent.PROGRESS,
                done=self.finished_jobs,
                total=len(self.workflow.dag),
            ),
        )


class JobRateLimiter:
    def __init__(self, limit: MaxJobsPerTimespan):
        self._limit: MaxJobsPerTimespan = limit
        self._jobs = deque()
        logger.debug(
            f"Submitting maximum {self._limit.max_jobs} job(s) over {self._limit.timespan} second(s)."
        )

    @property
    def max_jobs(self) -> int:
        return self._limit.max_jobs

    @property
    def timespan(self) -> int:
        return self._limit.timespan

    def register_jobs(self, n_jobs: int):
        currtime = time.time()
        self._jobs.extend(repeat(currtime, n_jobs))

    def get_free_jobs(self):
        # get the index of the last element that is older than the timespan
        index = bisect(self._jobs, time.time() - self.timespan)
        # remove the first index elements from the deque
        for _ in range(index):
            self._jobs.popleft()
        n_free = max(self.max_jobs - len(self._jobs), 0)
        return n_free
