__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os, signal, sys
import threading
import operator
from functools import partial
from collections import defaultdict
from itertools import chain, accumulate
from contextlib import ContextDecorator
import time

from snakemake.executors import DryrunExecutor, TouchExecutor, CPUExecutor
from snakemake.executors import (
    GenericClusterExecutor,
    SynchronousClusterExecutor,
    DRMAAExecutor,
    KubernetesExecutor,
    TibannaExecutor,
)
from snakemake.executors.google_lifesciences import GoogleLifeSciencesExecutor
from snakemake.exceptions import RuleException, WorkflowError, print_exception
from snakemake.shell import shell

from snakemake.logging import logger

from fractions import Fraction


def cumsum(iterable, zero=[0]):
    return list(chain(zero, accumulate(iterable)))


_ERROR_MSG_FINAL = (
    "Exiting because a job execution failed. " "Look above for error message"
)


class DummyRateLimiter(ContextDecorator):
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class JobScheduler:
    def __init__(
        self,
        workflow,
        dag,
        cores,
        local_cores=1,
        dryrun=False,
        touch=False,
        cluster=None,
        cluster_status=None,
        cluster_config=None,
        cluster_sync=None,
        drmaa=None,
        drmaa_log_dir=None,
        kubernetes=None,
        container_image=None,
        tibanna=None,
        tibanna_sfn=None,
        google_lifesciences=None,
        google_lifesciences_regions=None,
        google_lifesciences_location=None,
        google_lifesciences_cache=False,
        precommand="",
        tibanna_config=False,
        jobname=None,
        quiet=False,
        printreason=False,
        printshellcmds=False,
        keepgoing=False,
        max_jobs_per_second=None,
        max_status_checks_per_second=100,
        latency_wait=3,
        greediness=1.0,
        force_use_threads=False,
        assume_shared_fs=True,
        keepincomplete=False,
    ):
        """ Create a new instance of KnapsackJobScheduler. """
        from ratelimiter import RateLimiter

        self.cluster = cluster
        self.cluster_config = cluster_config
        self.cluster_sync = cluster_sync
        self.dag = dag
        self.workflow = workflow
        self.dryrun = dryrun
        self.touch = touch
        self.quiet = quiet
        self.keepgoing = keepgoing
        self.running = set()
        self.failed = set()
        self.finished_jobs = 0
        self.greediness = 1
        self.max_jobs_per_second = max_jobs_per_second
        self.keepincomplete = keepincomplete

        self.global_resources = {
            name: (sys.maxsize if res is None else res)
            for name, res in workflow.global_resources.items()
        }
        self.resources = dict(self.global_resources)

        use_threads = (
            force_use_threads
            or (os.name != "posix")
            or cluster
            or cluster_sync
            or drmaa
        )
        self._open_jobs = threading.Semaphore(0)
        self._lock = threading.Lock()

        self._errors = False
        self._finished = False
        self._job_queue = None
        self._submit_callback = self._noop
        self._finish_callback = partial(
            self._proceed,
            update_dynamic=not self.dryrun,
            print_progress=not self.quiet and not self.dryrun,
        )

        self._local_executor = None
        if dryrun:
            self._executor = DryrunExecutor(
                workflow,
                dag,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                latency_wait=latency_wait,
            )
        elif touch:
            self._executor = TouchExecutor(
                workflow,
                dag,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                latency_wait=latency_wait,
            )
        elif cluster or cluster_sync or (drmaa is not None):
            if not workflow.immediate_submit:
                # No local jobs when using immediate submit!
                # Otherwise, they will fail due to missing input
                self._local_executor = CPUExecutor(
                    workflow,
                    dag,
                    local_cores,
                    printreason=printreason,
                    quiet=quiet,
                    printshellcmds=printshellcmds,
                    latency_wait=latency_wait,
                    cores=local_cores,
                    keepincomplete=keepincomplete,
                )
            if cluster or cluster_sync:
                if cluster_sync:
                    constructor = SynchronousClusterExecutor
                else:
                    constructor = partial(
                        GenericClusterExecutor,
                        statuscmd=cluster_status,
                        max_status_checks_per_second=max_status_checks_per_second,
                    )

                self._executor = constructor(
                    workflow,
                    dag,
                    None,
                    submitcmd=(cluster or cluster_sync),
                    cluster_config=cluster_config,
                    jobname=jobname,
                    printreason=printreason,
                    quiet=quiet,
                    printshellcmds=printshellcmds,
                    latency_wait=latency_wait,
                    assume_shared_fs=assume_shared_fs,
                    keepincomplete=keepincomplete,
                )
                if workflow.immediate_submit:
                    self._submit_callback = partial(
                        self._proceed,
                        update_dynamic=False,
                        print_progress=False,
                        update_resources=False,
                        handle_job_success=False,
                    )
            else:
                self._executor = DRMAAExecutor(
                    workflow,
                    dag,
                    None,
                    drmaa_args=drmaa,
                    drmaa_log_dir=drmaa_log_dir,
                    jobname=jobname,
                    printreason=printreason,
                    quiet=quiet,
                    printshellcmds=printshellcmds,
                    latency_wait=latency_wait,
                    cluster_config=cluster_config,
                    assume_shared_fs=assume_shared_fs,
                    max_status_checks_per_second=max_status_checks_per_second,
                    keepincomplete=keepincomplete,
                )
        elif kubernetes:
            self._local_executor = CPUExecutor(
                workflow,
                dag,
                local_cores,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                latency_wait=latency_wait,
                cores=local_cores,
                keepincomplete=keepincomplete,
            )

            self._executor = KubernetesExecutor(
                workflow,
                dag,
                kubernetes,
                container_image=container_image,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                latency_wait=latency_wait,
                cluster_config=cluster_config,
                keepincomplete=keepincomplete,
            )
        elif tibanna:
            self._local_executor = CPUExecutor(
                workflow,
                dag,
                local_cores,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                use_threads=use_threads,
                latency_wait=latency_wait,
                cores=local_cores,
                keepincomplete=keepincomplete,
            )

            self._executor = TibannaExecutor(
                workflow,
                dag,
                cores,
                tibanna_sfn,
                precommand=precommand,
                tibanna_config=tibanna_config,
                container_image=container_image,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                latency_wait=latency_wait,
                keepincomplete=keepincomplete,
            )
        elif google_lifesciences:
            self._local_executor = CPUExecutor(
                workflow,
                dag,
                local_cores,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                latency_wait=latency_wait,
                cores=local_cores,
            )

            self._executor = GoogleLifeSciencesExecutor(
                workflow,
                dag,
                cores,
                container_image=container_image,
                regions=google_lifesciences_regions,
                location=google_lifesciences_location,
                cache=google_lifesciences_cache,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                latency_wait=latency_wait,
            )

        else:
            self._executor = CPUExecutor(
                workflow,
                dag,
                cores,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                use_threads=use_threads,
                latency_wait=latency_wait,
                cores=cores,
                keepincomplete=keepincomplete,
            )
        if self.max_jobs_per_second and not self.dryrun:
            max_jobs_frac = Fraction(self.max_jobs_per_second).limit_denominator()
            self.rate_limiter = RateLimiter(
                max_calls=max_jobs_frac.numerator, period=max_jobs_frac.denominator
            )

        else:
            # essentially no rate limit
            self.rate_limiter = DummyRateLimiter()

        self._user_kill = None
        signal.signal(signal.SIGTERM, self.exit_gracefully)
        self._open_jobs.release()

    @property
    def stats(self):
        try:
            return self._executor.stats
        except AttributeError:
            raise TypeError("Executor does not support stats")

    def candidate(self, job):
        """ Return whether a job is a candidate to be executed. """
        return (
            job not in self.running
            and job not in self.failed
            and (self.dryrun or (not job.dynamic_input and not self.dag.dynamic(job)))
        )

    @property
    def open_jobs(self):
        """ Return open jobs. """
        return filter(self.candidate, list(job for job in self.dag.ready_jobs))

    def schedule(self):
        """ Schedule jobs that are ready, maximizing cpu usage. """
        try:
            while True:
                # work around so that the wait does not prevent keyboard interrupts
                # while not self._open_jobs.acquire(False):
                #    time.sleep(1)
                self._open_jobs.acquire()

                # obtain needrun and running jobs in a thread-safe way
                with self._lock:
                    needrun = list(self.open_jobs)
                    running = list(self.running)
                    errors = self._errors
                    user_kill = self._user_kill

                # handle errors
                if user_kill or (not self.keepgoing and errors):
                    if user_kill == "graceful":
                        logger.info(
                            "Will exit after finishing " "currently running jobs."
                        )

                    if not running:
                        logger.info("Shutting down, this might take some time.")
                        self._executor.shutdown()
                        if not user_kill:
                            logger.error(_ERROR_MSG_FINAL)
                        return False
                    continue

                # normal shutdown because all jobs have been finished
                if not needrun and (not running or self.workflow.immediate_submit):
                    self._executor.shutdown()
                    if errors:
                        logger.error(_ERROR_MSG_FINAL)
                    return not errors

                # continue if no new job needs to be executed
                if not needrun:
                    continue

                # select jobs by solving knapsack problem (omit with dryrun)
                if self.dryrun:
                    run = needrun
                else:
                    logger.debug(
                        "Resources before job selection: {}".format(self.resources)
                    )
                    logger.debug(
                        "Ready jobs ({}):\n\t".format(len(needrun))
                        + "\n\t".join(map(str, needrun))
                    )
                    run = self.job_selector(needrun)
                    logger.debug(
                        "Selected jobs ({}):\n\t".format(len(run))
                        + "\n\t".join(map(str, run))
                    )
                    logger.debug(
                        "Resources after job selection: {}".format(self.resources)
                    )
                # update running jobs
                with self._lock:
                    self.running.update(run)

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
        for name, value in job.resources.items():
            if name in self.resources:
                value = self.calc_resource(name, value)
                self.resources[name] += value

    def _proceed(
        self,
        job,
        update_dynamic=True,
        print_progress=False,
        update_resources=True,
        handle_job_success=True,
    ):
        """ Do stuff after job is finished. """
        with self._lock:
            if handle_job_success:
                # by calling this behind the lock, we avoid race conditions
                try:
                    self.get_executor(job).handle_job_success(job)
                except (RuleException, WorkflowError) as e:
                    # if an error occurs while processing job output,
                    # we do the same as in case of errors during execution
                    print_exception(e, self.workflow.linemaps)
                    self._handle_error(job)
                    return

            try:
                self.dag.finish(job, update_dynamic=update_dynamic)
            except (RuleException, WorkflowError) as e:
                # if an error occurs while processing job output,
                # we do the same as in case of errors during execution
                print_exception(e, self.workflow.linemaps)
                self._handle_error(job)
                return

            if update_resources:
                # normal jobs have len=1, group jobs have len>1
                self.finished_jobs += len(job)
                self.running.remove(job)
                self._free_resources(job)

            if print_progress:
                if job.is_group():
                    for j in job:
                        logger.job_finished(jobid=j.jobid)
                else:
                    logger.job_finished(jobid=job.jobid)
                self.progress()

            if (
                any(self.open_jobs)
                or not self.running
                or self.workflow.immediate_submit
            ):
                # go on scheduling if open jobs are ready or no job is running
                self._open_jobs.release()

    def _error(self, job):
        with self._lock:
            self._handle_error(job)

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
            logger.info("Trying to restart job {}.".format(self.dag.jobid(job)))
            job.attempt += 1
        else:
            self._errors = True
            self.failed.add(job)
            if self.keepgoing:
                logger.info("Job failed, going on with independent jobs.")
        self._open_jobs.release()

    def exit_gracefully(self, *args):
        with self._lock:
            self._user_kill = "graceful"
        self._open_jobs.release()

    def job_selector(self, jobs):
        """
        Using the greedy heuristic from
        "A Greedy Algorithm for the General Multidimensional Knapsack
Problem", Akcay, Li, Xu, Annals of Operations Research, 2012

        Args:
            jobs (list):    list of jobs
        """
        with self._lock:
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

            solution = [job for job, sel in zip(jobs, x) if sel]
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
                "jobs are connected via a pipe and have to run "
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
        res = job.resources
        return [
            self.calc_resource(name, res.get(name, 0)) for name in self.global_resources
        ]

    def job_reward(self, job):
        if self.touch or self.dryrun or self.workflow.immediate_submit:
            temp_size = 0
            input_size = 0
        else:
            temp_size = self.dag.temp_size(job)
            input_size = job.inputsize

        # Usually, this should guide the scheduler to first schedule all jobs
        # that remove the largest temp file, then the second largest and so on.
        # Since the weight is summed up, it can in theory be that it sometimes
        # prefers a set of many jobs that all depend on smaller temp files though.
        # A real solution to the problem is therefore to use dummy jobs that
        # ensure selection of groups of jobs that together delete the same temp
        # file.

        return (job.priority, temp_size, input_size)

    def progress(self):
        """ Display the progress. """
        logger.progress(done=self.finished_jobs, total=len(self.dag))
