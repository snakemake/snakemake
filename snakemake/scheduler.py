__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os, signal
import threading
import operator
from functools import partial
from collections import defaultdict
from itertools import chain, accumulate

from snakemake.executors import DryrunExecutor, TouchExecutor, CPUExecutor
from snakemake.executors import GenericClusterExecutor, SynchronousClusterExecutor, DRMAAExecutor

from snakemake.logging import logger


def cumsum(iterable, zero=[0]):
    return list(chain(zero, accumulate(iterable)))


_ERROR_MSG_FINAL = ("Exiting because a job execution failed. "
                    "Look above for error message")


class JobScheduler:
    def __init__(self, workflow, dag, cores,
                 local_cores=1,
                 dryrun=False,
                 touch=False,
                 cluster=None,
                 cluster_config=None,
                 cluster_sync=None,
                 drmaa=None,
                 drmaa_log_dir=None,
                 jobname=None,
                 quiet=False,
                 printreason=False,
                 printshellcmds=False,
                 keepgoing=False,
                 max_jobs_per_second=None,
                 latency_wait=3,
                 benchmark_repeats=1,
                 greediness=1.0,
                 force_use_threads=False):
        """ Create a new instance of KnapsackJobScheduler. """
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

        self.resources = dict(self.workflow.global_resources)

        use_threads = force_use_threads or (os.name != "posix") or cluster or cluster_sync or drmaa
        self._open_jobs = threading.Event()
        self._lock = threading.Lock()

        self._errors = False
        self._finished = False
        self._job_queue = None
        self._submit_callback = self._noop
        self._finish_callback = partial(
            self._proceed,
            update_dynamic=not self.dryrun,
            print_progress=not self.quiet and not self.dryrun)

        self._local_executor = None
        if dryrun:
            self._executor = DryrunExecutor(workflow, dag,
                                            printreason=printreason,
                                            quiet=quiet,
                                            printshellcmds=printshellcmds,
                                            latency_wait=latency_wait)
            self.job_reward = self.dryrun_job_reward
        elif touch:
            self._executor = TouchExecutor(workflow, dag,
                                           printreason=printreason,
                                           quiet=quiet,
                                           printshellcmds=printshellcmds,
                                           latency_wait=latency_wait)
        elif cluster or cluster_sync or (drmaa is not None):
            workers = min(max(1, sum(1 for _ in dag.local_needrun_jobs)), local_cores)
            self._local_executor = CPUExecutor(
                workflow, dag, workers,
                printreason=printreason,
                quiet=quiet,
                printshellcmds=printshellcmds,
                use_threads=use_threads,
                latency_wait=latency_wait,
                benchmark_repeats=benchmark_repeats,
                cores=local_cores)
            if cluster or cluster_sync:
                constructor = SynchronousClusterExecutor if cluster_sync \
                              else GenericClusterExecutor
                self._executor = constructor(
                    workflow, dag, None,
                    submitcmd=(cluster or cluster_sync),
                    cluster_config=cluster_config,
                    jobname=jobname,
                    printreason=printreason,
                    quiet=quiet,
                    printshellcmds=printshellcmds,
                    latency_wait=latency_wait,
                    benchmark_repeats=benchmark_repeats,
                    max_jobs_per_second=max_jobs_per_second)
                if workflow.immediate_submit:
                    self.job_reward = self.dryrun_job_reward
                    self._submit_callback = partial(self._proceed,
                                                    update_dynamic=False,
                                                    print_progress=False,
                                                    update_resources=False, )
            else:
                self._executor = DRMAAExecutor(
                    workflow, dag, None,
                    drmaa_args=drmaa,
                    drmaa_log_dir=drmaa_log_dir,
                    jobname=jobname,
                    printreason=printreason,
                    quiet=quiet,
                    printshellcmds=printshellcmds,
                    latency_wait=latency_wait,
                    benchmark_repeats=benchmark_repeats,
                    cluster_config=cluster_config,
                    max_jobs_per_second=max_jobs_per_second)
        else:
            # local execution or execution of cluster job
            # calculate how many parallel workers the executor shall spawn
            # each job has at least one thread, hence we need to have
            # the minimum of given cores and number of jobs
            workers = min(cores, max(1, len(dag)))
            self._executor = CPUExecutor(workflow, dag, workers,
                                         printreason=printreason,
                                         quiet=quiet,
                                         printshellcmds=printshellcmds,
                                         use_threads=use_threads,
                                         latency_wait=latency_wait,
                                         benchmark_repeats=benchmark_repeats,
                                         cores=cores)
        self._open_jobs.set()

    @property
    def stats(self):
        try:
            return self._executor.stats
        except AttributeError:
            raise TypeError("Executor does not support stats")

    def candidate(self, job):
        """ Return whether a job is a candidate to be executed. """
        return (job not in self.running and job not in self.failed and
                (self.dryrun or
                 (not job.dynamic_input and not self.dag.dynamic(job))))

    @property
    def open_jobs(self):
        """ Return open jobs. """
        return filter(self.candidate, list(self.dag.ready_jobs))

    def schedule(self):
        """ Schedule jobs that are ready, maximizing cpu usage. """
        try:
            while True:
                # work around so that the wait does not prevent keyboard interrupts
                while not self._open_jobs.wait(1):
                    pass

                # obtain needrun and running jobs in a thread-safe way
                with self._lock:
                    needrun = list(self.open_jobs)
                    running = list(self.running)
                # free the event
                self._open_jobs.clear()

                # handle errors
                if not self.keepgoing and self._errors:
                    logger.info("Will exit after finishing "
                                "currently running jobs.")
                    if not running:
                        self._executor.shutdown()
                        logger.error(_ERROR_MSG_FINAL)
                        return False
                    continue
                # normal shutdown because all jobs have been finished
                if not needrun and not running:
                    self._executor.shutdown()
                    if self._errors:
                        logger.error(_ERROR_MSG_FINAL)
                    return not self._errors

                # continue if no new job needs to be executed
                if not needrun:
                    continue

                logger.debug("Resources before job selection: {}".format(
                    self.resources))
                logger.debug("Ready jobs ({}):\n\t".format(len(needrun)) +
                             "\n\t".join(map(str, needrun)))

                # select jobs by solving knapsack problem
                run = self.job_selector(needrun)
                logger.debug("Selected jobs ({}):\n\t".format(len(run)) +
                             "\n\t".join(map(str, run)))
                # update running jobs
                with self._lock:
                    self.running.update(run)
                logger.debug(
                    "Resources after job selection: {}".format(self.resources))
                # actually run jobs
                for job in run:
                    self.run(job)
        except (KeyboardInterrupt, SystemExit):
            logger.info("Terminating processes on user request.")
            self._executor.cancel()
            with self._lock:
                running = list(self.running)
            for job in running:
                job.cleanup()
            return False

    def get_executor(self, job):
        if self._local_executor is None:
            return self._executor
        else:
            return self._local_executor if self.workflow.is_local(
                job.rule) else self._executor

    def run(self, job):
        self.get_executor(job).run(job,
            callback=self._finish_callback,
            submit_callback=self._submit_callback,
            error_callback=self._error)

    def _noop(self, job):
        pass

    def _free_resources(self, job):
        for name, value in job.resources.items():
            if name in self.resources:
                value = self.calc_resource(name, value)
                self.resources[name] += value
                logger.debug("Releasing {} {} (now {}).".format(
                    value, name, self.resources[name]))

    def _proceed(self, job,
                 update_dynamic=True,
                 print_progress=False,
                 update_resources=True):
        """ Do stuff after job is finished. """
        with self._lock:
            # by calling this behind the lock, we avoid race conditions
            self.get_executor(job).handle_job_success(job)
            self.dag.finish(job, update_dynamic=update_dynamic)

            if update_resources:
                self.finished_jobs += 1
                self.running.remove(job)
                self._free_resources(job)


            if print_progress:
                logger.job_finished(jobid=self.dag.jobid(job))
                self.progress()

            if any(self.open_jobs) or not self.running:
                # go on scheduling if open jobs are ready or no job is running
                self._open_jobs.set()

    def _error(self, job):
        """Clear jobs and stop the workflow.

        If Snakemake is configured to restart jobs then the job might have
        "restart_times" left and we just decrement and let the scheduler
        try to run the job again.
        """
        with self._lock:
            self.get_executor(job).handle_job_error(job)
            self.running.remove(job)
            self._free_resources(job)
            self._open_jobs.set()
            if job.restart_times > 0:
                msg = (
                    ("Trying to restart job for rule {} with "
                     "wildcards {}").format(
                         job.rule.name, job.wildcards_dict))
                logger.info(msg
                    )
                job.restart_times -= 1
            else:
                self._errors = True
                self.failed.add(job)
                if self.keepgoing:
                    logger.info("Job failed, going on with independent jobs.")

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

            b = [self.resources[name]
                 for name in self.workflow.global_resources
                 ]  # resource capacities

            while True:
                # Step 2: compute effective capacities
                y = [(min((min(u[j], b_i // a_j_i) if a_j_i > 0 else u[j])
                          for b_i, a_j_i in zip(b, a[j]) if a_j_i) if j in E
                      else 0) for j in range(n)]
                if not any(y):
                    break
                y = [(max(1, int(self.greediness * y_j)) if y_j > 0 else 0)
                     for y_j in y]

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
            for name, b_i in zip(self.workflow.global_resources, b):
                self.resources[name] = b_i
            return solution

    def calc_resource(self, name, value):
        return min(value, self.workflow.global_resources[name])

    def rule_weight(self, rule):
        res = rule.resources
        return [self.calc_resource(name, res.get(name, 0))
                for name in self.workflow.global_resources]

    def job_weight(self, job):
        res = job.resources
        return [self.calc_resource(name, res.get(name, 0))
                for name in self.workflow.global_resources]

    def job_reward(self, job):
        return (self.dag.priority(job), self.dag.temp_input_count(job), self.dag.downstream_size(job),
                0 if self.touch else job.inputsize)

    def dryrun_job_reward(self, job):
        return (self.dag.priority(job), self.dag.temp_input_count(job), self.dag.downstream_size(job))

    def progress(self):
        """ Display the progress. """
        logger.progress(done=self.finished_jobs, total=len(self.dag))
