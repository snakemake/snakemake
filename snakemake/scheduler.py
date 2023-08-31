__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os, signal, sys
import threading

from functools import partial
from itertools import chain, accumulate
from contextlib import ContextDecorator

from snakemake_interface_executor_plugins.scheduler import JobSchedulerExecutorInterface
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
from snakemake_interface_executor_plugins.registry import Plugin as ExecutorPlugin

from snakemake.exceptions import RuleException, WorkflowError, print_exception
from snakemake.logging import logger

from fractions import Fraction

registry = ExecutorPluginRegistry()


def cumsum(iterable, zero=[0]):
    return list(chain(zero, accumulate(iterable)))


_ERROR_MSG_FINAL = (
    "Exiting because a job execution failed. Look above for error message"
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
    def __init__(
        self,
        workflow,
        executor_plugin: ExecutorPlugin
    ):
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
        self.max_jobs_per_second = self.workflow.scheduling_settings.max_jobs_per_second
        self._tofinish = []
        self._toerror = []
        self.handle_job_success = True
        self.update_resources = True
        self.print_progress = not self.quiet and not self.dryrun
        self.update_dynamic = not self.dryrun

        self.global_resources = {
            name: (sys.maxsize if res is None else res)
            for name, res in workflow.global_resources.items()
        }

        if workflow.global_resources["_nodes"] is not None:
            # Do not restrict cores locally if nodes are used (i.e. in case of cluster/cloud submission).
            self.global_resources["_cores"] = sys.maxsize
        self.resources = dict(self.global_resources)

        self._open_jobs = threading.Semaphore(0)
        self._lock = threading.Lock()

        self._errors = False
        self._executor_error = None
        self._finished = False
        self._job_queue = None
        self._last_job_selection_empty = False
        self._submit_callback = self._noop
        self._finish_callback = self._proceed

        self._local_executor = None

        if self.workflow.executor_plugin.common_settings.local_exec:
            self._executor = executor_plugin.executor(
                self.workflow,
                logger,
            )
        else:
            self._executor = executor_plugin.executor(
                self.workflow,
                logger,
            )
            self._local_executor = ExecutorPluginRegistry().get_executor("local").executor(
                self.workflow,
                logger,
            )

        # elif slurm:
        #     if ON_WINDOWS:
        #         raise WorkflowError("SLURM execution is not supported on Windows.")
        #     self._local_executor = CPUExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         local_cores,
        #     )
        #     # we need to adjust the maximum status checks per second
        #     # on a SLURM cluster, to not overstrain the scheduler;
        #     # timings for tested SLURM clusters, extracted from --verbose
        #     # output with:
        #     # ```
        #     #   grep "sacct output" .snakemake/log/2023-02-13T210004.601290.snakemake.log | \
        #     #   awk '{ counter += 1; sum += $6; sum_of_squares += ($6)^2 } \
        #     #     END { print "average: ",sum/counter," sd: ",sqrt((sum_of_squares - sum^2/counter) / counter); }
        #     # ````
        #     #   * cluster 1:
        #     #     * sacct:    average:  0.073896   sd:  0.0640178
        #     #     * scontrol: average:  0.0193017  sd:  0.0358858
        #     # Thus, 2 status checks per second should leave enough
        #     # capacity for everybody.
        #     # TODO: check timings on other slurm clusters, to:
        #     #   * confirm that this cap is reasonable
        #     #   * check if scontrol is the quicker option across the board
        #     if max_status_checks_per_second > 2:
        #         max_status_checks_per_second = 2

        #     self._executor = SlurmExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         max_status_checks_per_second=max_status_checks_per_second,
        #     )

        # elif slurm_jobstep:
        #     self._executor = SlurmJobstepExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #     )
        #     self._local_executor = self._executor

        # elif cluster or cluster_sync or (drmaa is not None):
        #     if not workflow.remote_execution_settings.immediate_submit:
        #         # No local jobs when using immediate submit!
        #         # Otherwise, they will fail due to missing input
        #         self._local_executor = CPUExecutor(
        #             workflow,
        #             dag,
        #             self.stats,
        #             logger,
        #             local_cores,
        #         )

        #     if cluster or cluster_sync:
        #         if cluster_sync:
        #             constructor = SynchronousClusterExecutor
        #         else:
        #             constructor = partial(
        #                 GenericClusterExecutor,
        #                 statuscmd=cluster_status,
        #                 cancelcmd=cluster_cancel,
        #                 cancelnargs=cluster_cancel_nargs,
        #                 sidecarcmd=cluster_sidecar,
        #                 max_status_checks_per_second=max_status_checks_per_second,
        #             )

        #         self._executor = constructor(
        #             workflow,
        #             dag,
        #             self.stats,
        #             logger,
        #             submitcmd=(cluster or cluster_sync),
        #             jobname=jobname,
        #         )
        #         if workflow.remote_execution_settings.immediate_submit:
        #             self._submit_callback = self._proceed
        #             self.update_dynamic = False
        #             self.print_progress = False
        #             self.update_resources = False
        #             self.handle_job_success = False
        #     else:
        #         self._executor = DRMAAExecutor(
        #             workflow,
        #             dag,
        #             self.stats,
        #             logger,
        #             drmaa_args=drmaa,
        #             drmaa_log_dir=drmaa_log_dir,
        #             jobname=jobname,
        #             max_status_checks_per_second=max_status_checks_per_second,
        #         )
        # elif kubernetes:
        #     self._local_executor = CPUExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         local_cores,
        #     )

        #     self._executor = KubernetesExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         kubernetes,
        #         container_image=container_image,
        #         k8s_cpu_scalar=k8s_cpu_scalar,
        #         k8s_service_account_name=k8s_service_account_name,
        #     )
        # elif tibanna:
        #     self._local_executor = CPUExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         local_cores,
        #         use_threads=use_threads,
        #     )

        #     self._executor = TibannaExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         cores,
        #         tibanna_sfn,
        #         precommand=precommand,
        #         tibanna_config=tibanna_config,
        #         container_image=container_image,
        #     )

        # elif flux:
        #     self._local_executor = CPUExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         local_cores,
        #     )

        #     self._executor = FluxExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #     )

        # elif az_batch:
        #     try:
        #         from snakemake.executors.azure_batch import AzBatchExecutor
        #     except ImportError as e:
        #         raise WorkflowError(
        #             "Unable to load Azure Batch executor. You have to install "
        #             "the msrest, azure-core, azure-batch, azure-mgmt-batch, and azure-identity packages.",
        #             e,
        #         )
        #     self._local_executor = CPUExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         local_cores,
        #     )
        #     self._executor = AzBatchExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         container_image=container_image,
        #         az_batch_account_url=az_batch_account_url,
        #         az_batch_enable_autoscale=az_batch_enable_autoscale,
        #     )

        # elif google_lifesciences:
        #     self._local_executor = CPUExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         local_cores,
        #     )

        #     self._executor = GoogleLifeSciencesExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         container_image=container_image,
        #         regions=google_lifesciences_regions,
        #         location=google_lifesciences_location,
        #         cache=google_lifesciences_cache,
        #         service_account_email=google_lifesciences_service_account_email,
        #         network=google_lifesciences_network,
        #         subnetwork=google_lifesciences_subnetwork,
        #         preemption_default=preemption_default,
        #         preemptible_rules=preemptible_rules,
        #     )
        # elif tes:
        #     self._local_executor = CPUExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         local_cores,
        #     )

        #     self._executor = TaskExecutionServiceExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         tes_url=tes,
        #         container_image=container_image,
        #     )

        # else:
        #     self._executor = CPUExecutor(
        #         workflow,
        #         dag,
        #         self.stats,
        #         logger,
        #         cores,
        #         use_threads=use_threads,
        #     )
        from throttler import Throttler

        if not self.dryrun:
            max_jobs_frac = Fraction(self.max_jobs_per_second).limit_denominator()
            self.rate_limiter = Throttler(
                rate_limit=max_jobs_frac.numerator, period=max_jobs_frac.denominator
            )
        else:
            # essentially no rate limit
            self.rate_limiter = DummyRateLimiter()

        # Choose job selector (greedy or ILP)
        self.job_selector = self.job_selector_greedy
        if self.workflow.scheduling_settings.scheduler == "ilp":
            import pulp

            if pulp.apis.LpSolverDefault is None:
                logger.warning(
                    "Falling back to greedy scheduler because no default "
                    "solver is found for pulp (you have to install either "
                    "coincbc or glpk)."
                )
            else:
                self.job_selector = self.job_selector_ilp

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
        jobs = self.workflow.dag.ready_jobs

        if not self.dryrun:
            jobs = [
                job
                for job in jobs
                if not job.dynamic_input and not self.workflow.dag.dynamic(job)
            ]
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
                        return False
                    continue

                # all runnable jobs have finished, normal shutdown
                if not needrun and (not running or self.workflow.remote_execution_settings.immediate_submit):
                    self._executor.shutdown()
                    if errors:
                        logger.error(_ERROR_MSG_FINAL)
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
                    logger.debug(
                        f"Ready jobs ({len(needrun)})"
                        # + "\n\t".join(map(str, needrun))
                    )

                    if not self._last_job_selection_empty:
                        logger.info("Select jobs to execute...")
                    run = self.job_selector(needrun)
                    self._last_job_selection_empty = not run

                    logger.debug(
                        f"Selected jobs ({len(run)})"
                        # + "\n\t".join(map(str, run))
                    )
                    logger.debug(f"Resources after job selection: {self.resources}")

                # update running jobs
                with self._lock:
                    self.running.update(run)
                    # remove from ready_jobs
                    self.workflow.dag.register_running(run)

                # actually run jobs
                local_runjobs = [job for job in run if job.is_local]
                runjobs = [job for job in run if not job.is_local]
                self.run(local_runjobs, executor=self._local_executor or self._executor)
                self.run(runjobs)
        except (KeyboardInterrupt, SystemExit):
            logger.info(
                "Terminating processes on user request, this might take some time."
            )
            self._executor.cancel()
            return False

    def _finish_jobs(self):
        # must be called from within lock
        # clear the global tofinish such that parallel calls do not interfere
        for job in self._tofinish:
            if self.handle_job_success:
                try:
                    self.get_executor(job).handle_job_success(job)
                except (RuleException, WorkflowError) as e:
                    # if an error occurs while processing job output,
                    # we do the same as in case of errors during execution
                    print_exception(e, self.workflow.linemaps)
                    self._handle_error(job)
                    continue

            if self.update_resources:
                # normal jobs have len=1, group jobs have len>1
                self.finished_jobs += len(job)
                self.running.remove(job)
                self._free_resources(job)

            if self.print_progress:
                if job.is_group():
                    for j in job:
                        logger.job_finished(jobid=j.jobid)
                else:
                    logger.job_finished(jobid=job.jobid)
                self.progress()

            self.workflow.dag.finish(job, update_dynamic=self.update_dynamic)
        self._tofinish.clear()

    def _error_jobs(self):
        # must be called from within lock
        for job in self._toerror:
            self._handle_error(job)
        self._toerror.clear()

    def run(self, jobs, executor=None):
        if executor is None:
            executor = self._executor

        executor.run_jobs(
            jobs,
            callback=self._finish_callback,
            submit_callback=self._submit_callback,
            error_callback=self._error,
        )

    def get_executor(self, job):
        if job.is_local and self._local_executor is not None:
            return self._local_executor
        return self._executor

    def _noop(self, job):
        pass

    def _free_resources(self, job):
        for name, value in job.scheduler_resources.items():
            if name in self.resources:
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

    def _error(self, job):
        with self._lock:
            self._toerror.append(job)
            self._open_jobs.release()

    def _handle_error(self, job):
        """Clear jobs and stop the workflow.

        If Snakemake is configured to restart jobs then the job might have
        "restart_times" left and we just decrement and let the scheduler
        try to run the job again.
        """
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

    def job_selector_ilp(self, jobs):
        """
        Job scheduling by optimization of resource usage by solving ILP using pulp
        """
        import pulp
        from pulp import lpSum
        from stopit import ThreadingTimeout as Timeout, TimeoutException

        if len(jobs) == 1:
            logger.debug(
                "Using greedy selector because only single job has to be scheduled."
            )
            return self.job_selector_greedy(jobs)

        with self._lock:
            if not self.resources["_cores"]:
                return set()

            # assert self.resources["_cores"] > 0
            scheduled_jobs = {
                job: pulp.LpVariable(
                    f"job_{idx}", lowBound=0, upBound=1, cat=pulp.LpInteger
                )
                for idx, job in enumerate(jobs)
            }

            def size_gb(f):
                if self.touch:
                    # In case of touch mode, there is no need to prioritize based on size.
                    # We cannot access it anyway, because the files might be temporary and
                    # not present.
                    return 0
                else:
                    return f.size / 1e9

            temp_files = {
                temp_file for job in jobs for temp_file in self.workflow.dag.temp_input(job)
            }

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
                sum([size_gb(temp_file) for temp_file in temp_files]), 1
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
                        temp_file_deletable[temp_file] * size_gb(temp_file)
                        for temp_file in temp_files
                    ]
                )
                + lpSum(
                    [
                        temp_job_improvement[temp_file] * size_gb(temp_file)
                        for temp_file in temp_files
                    ]
                )
            )

            # Constraints:
            for name in self.workflow.global_resources:
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

        try:
            with Timeout(10, swallow_exc=False):
                self._solve_ilp(prob)
        except TimeoutException as e:
            logger.warning(
                "Failed to solve scheduling problem with ILP solver in time (10s). "
                "Falling back to greedy solver."
            )
            return self.job_selector_greedy(jobs)
        except pulp.apis.core.PulpSolverError as e:
            logger.warning(
                "Failed to solve scheduling problem with ILP solver. Falling back to greedy solver. "
                "Run Snakemake with --verbose to see the full solver output for debugging the problem."
            )
            return self.job_selector_greedy(jobs)

        selected_jobs = set(
            job for job, variable in scheduled_jobs.items() if variable.value() == 1.0
        )

        if not selected_jobs:
            # No selected jobs. This could be due to insufficient resources or a failure in the ILP solver
            # Hence, we silently fall back to the greedy solver to make sure that we don't miss anything.
            return self.job_selector_greedy(jobs)

        for name in self.workflow.global_resources:
            self.resources[name] -= sum(
                [job.scheduler_resources.get(name, 0) for job in selected_jobs]
            )
        return selected_jobs

    def _solve_ilp(self, prob):
        import pulp

        old_path = os.environ["PATH"]
        if self.workflow.scheduling_settings.scheduler_solver_path is None:
            # Temporarily prepend the given snakemake env to the path, such that the solver can be found in any case.
            # This is needed for cluster envs, where the cluster job might have a different environment but
            # still needs access to the solver binary.
            os.environ["PATH"] = "{}:{}".format(
                self.workflow.scheduling_settings.scheduler_solver_path, os.environ["PATH"]
            )
        try:
            solver = (
                pulp.get_solver(self.workflow.scheduling_settings.scheduler_ilp_solver)
                if self.workflow.scheduling_settings.scheduler_ilp_solver
                else pulp.apis.LpSolverDefault
            )
        finally:
            os.environ["PATH"] = old_path
        solver.msg = self.workflow.output_settings.verbose
        prob.solve(solver)

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
        with self._lock:
            if not self.resources["_cores"]:
                return set()
            # each job is an item with one copy (0-1 MDKP)
            n = len(jobs)
            x = [0] * n  # selected jobs
            E = set(range(n))  # jobs still free to select
            u = [1] * n
            a = list(map(self.job_weight, jobs))  # resource usage of jobs
            c = list(map(self.job_reward, jobs))  # job rewards

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
            # update resources
            for name, b_i in zip(self.global_resources, b):
                self.resources[name] = b_i
            return solution

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

    def job_reward(self, job):
        if self.touch or self.dryrun or self.workflow.remote_execution_settings.immediate_submit:
            temp_size = 0
            input_size = 0
        else:
            try:
                temp_size = self.workflow.dag.temp_size(job)
                input_size = job.inputsize
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
        logger.progress(done=self.finished_jobs, total=len(self.workflow.dag))
