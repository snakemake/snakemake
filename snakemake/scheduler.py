# -*- coding: utf-8 -*-

import os, signal
import threading
import multiprocessing
import operator
from functools import partial
from collections import defaultdict
from itertools import chain, accumulate

from snakemake.executors import DryrunExecutor, TouchExecutor
from snakemake.executors import GenericClusterExecutor, CPUExecutor, DRMAAExecutor
from snakemake.logging import logger

__author__ = "Johannes Köster"



def cumsum(iterable, zero=[0]):
    return list(chain(zero, accumulate(iterable)))

_ERROR_MSG_FINAL = (
    "Exiting because a job execution failed. "
    "Look above for error message")


class JobScheduler:
    def __init__(
        self,
        workflow,
        dag,
        cores,
        dryrun=False,
        touch=False,
        cluster=None,
        drmaa=None,
        jobname=None,
        immediate_submit=False,
        quiet=False,
        printreason=False,
        printshellcmds=False,
        keepgoing=False,
        latency_wait=3,
        benchmark_repeats=3,
        greedyness=1.0
    ):
        """ Create a new instance of KnapsackJobScheduler. """
        self.cluster = cluster
        self.dag = dag
        self.workflow = workflow
        self.dryrun = dryrun
        self.quiet = quiet
        self.keepgoing = keepgoing
        self.running = set()
        self.failed = set()
        self.finished_jobs = 0
        self.greedyness = greedyness
        self.select_by_rule = False
        if not self.select_by_rule:
            self.greedyness = 1

        self.resources = dict(self.workflow.global_resources)

        use_threads = os.name != "posix"
        if not use_threads:
            self._open_jobs = multiprocessing.Event()
            self._lock = multiprocessing.Lock()
        else:
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

        if dryrun:
            self._executor = DryrunExecutor(
                workflow, dag, printreason=printreason,
                quiet=quiet, printshellcmds=printshellcmds,
                latency_wait=latency_wait)
            self.rule_reward = self.dryrun_rule_reward
            self.job_reward = self.dryrun_job_reward
        elif touch:
            self._executor = TouchExecutor(
                workflow, dag, printreason=printreason,
                quiet=quiet, printshellcmds=printshellcmds,
                latency_wait=latency_wait)
        elif cluster or (drmaa is not None):
            # TODO properly set cores
            self._local_executor = CPUExecutor(
                workflow, dag, multiprocessing.cpu_count(), printreason=printreason,
                quiet=quiet, printshellcmds=printshellcmds,
                threads=use_threads,
                latency_wait=latency_wait, benchmark_repeats=benchmark_repeats)
            if cluster:
                self._executor = GenericClusterExecutor(
                    workflow, dag, None, submitcmd=cluster, jobname=jobname,
                    printreason=printreason, quiet=quiet,
                    printshellcmds=printshellcmds, latency_wait=latency_wait,
                    benchmark_repeats=benchmark_repeats)
                if immediate_submit:
                    self.rule_reward = self.dryrun_rule_reward
                    self.job_reward = self.dryrun_job_reward
                    self._submit_callback = partial(
                        self._proceed,
                        update_dynamic=False,
                        print_progress=False,
                        update_resources=False)
                else:
                    self.run = self.run_cluster_or_local
            else:
                self._executor = DRMAAExecutor(
                    workflow, dag, None, drmaa_args=drmaa, jobname=jobname,
                    printreason=printreason, quiet=quiet,
                    printshellcmds=printshellcmds, latency_wait=latency_wait,
                    benchmark_repeats=benchmark_repeats)
                self.run = self.run_cluster_or_local
        else:
            self._executor = CPUExecutor(
                workflow, dag, cores, printreason=printreason,
                quiet=quiet, printshellcmds=printshellcmds,
                threads=use_threads,
                latency_wait=latency_wait,
                benchmark_repeats=benchmark_repeats)
        self._open_jobs.set()

    @property
    def stats(self):
        try:
            return self._executor.stats
        except AttributeError:
            raise TypeError("Executor does not support stats")

    def candidate(self, job):
        """ Return whether a job is a candidate to be executed. """
        return (job not in self.running and job not in self.failed and (self.dryrun or
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

                self._open_jobs.clear()
                if not self.keepgoing and self._errors:
                    logger.info("Will exit after finishing "
                        "currently running jobs.")
                    if not self.running:
                        self._executor.shutdown()
                        logger.error(_ERROR_MSG_FINAL)
                        return False
                    continue
                if not any(self.open_jobs) and not self.running:
                    self._executor.shutdown()
                    if self._errors:
                        logger.error(_ERROR_MSG_FINAL)
                    return not self._errors

                needrun = list(self.open_jobs)
                if not needrun:
                    continue

                logger.debug("Resources before job selection: {}".format(self.resources))
                logger.debug("Ready jobs ({}):\n\t".format(len(needrun)) + "\n\t".join(map(str, needrun)))

                run = self.job_selector(needrun)
                logger.debug("Selected jobs ({}):\n\t".format(len(run)) + "\n\t".join(map(str, run)))
                self.running.update(run)
                logger.debug("Resources after job selection: {}".format(self.resources))
                for job in run:
                    self.run(job)
        except (KeyboardInterrupt, SystemExit):
            logger.info("Terminating processes on user request.")
            self._executor.cancel()
            for job in self.running:
                job.cleanup()
            return False

    def run(self, job):
        self._executor.run(
            job, callback=self._finish_callback,
            submit_callback=self._submit_callback,
            error_callback=self._error)

    def run_cluster_or_local(self, job):
        executor = self._local_executor if self.workflow.is_local(job.rule) else self._executor
        executor.run(
            job, callback=self._finish_callback,
            submit_callback=self._submit_callback,
            error_callback=self._error)

    def _noop(self, job):
        pass

    def _free_resources(self, job):
        for name, value in job.resources.items():
            if name in self.resources:
                value = self.calc_resource(name, value)
                self.resources[name] += value
                logger.debug("Releasing {} {} (now {}).".format(value, name, self.resources[name]))

    def _proceed(
        self, job, update_dynamic=True, print_progress=False,
        update_resources=True):
        """ Do stuff after job is finished. """
        with self._lock:
            if update_resources:
                self.finished_jobs += 1
                self.running.remove(job)
                self._free_resources(job)

            self.dag.finish(job, update_dynamic=update_dynamic)

            logger.job_finished(jobid=self.dag.jobid(job))

            if print_progress:
                self.progress()

            if any(self.open_jobs) or not self.running:
                # go on scheduling if open jobs are ready or no job is running
                self._open_jobs.set()

    def _error(self, job):
        """ Clear jobs and stop the workflow. """
        with self._lock:
            self._errors = True
            self.running.remove(job)
            self.failed.add(job)
            self._free_resources(job)
            if self.keepgoing:
                logger.info("Job failed, going on with independent jobs.")
            self._open_jobs.set()

    def job_selector(self, jobs):
        """
        Using the greedy heuristic from
        "A Greedy Algorithm for the General Multidimensional Knapsack
Problem", Akcay, Li, Xu, Annals of Operations Research, 2012

        Args:
            jobs (list):    list of jobs
        """
        with self._lock:
            if self.select_by_rule:
                # solve over the rules instead of jobs (much less, but might miss the best solution)
                # each rule is an item with as many copies as jobs
                _jobs = defaultdict(list)
                for job in jobs:
                    _jobs[job.rule].append(job)

                jobs = _jobs

                # sort the jobs by priority
                for _jobs in jobs.values():
                    _jobs.sort(key=self.dag.priority, reverse=True)
                rules = list(jobs)

                # Step 1: initialization
                n = len(rules)
                x = [0] * n  # selected jobs of each rule
                E = set(range(n))  # rules free to select
                u = [len(jobs[rule]) for rule in rules]  # number of jobs left
                a = list(map(self.rule_weight, rules))  # resource usage of rules
                c = list(map(partial(self.rule_reward, jobs=jobs), rules))  # matrix of cumulative rewards over jobs

                def calc_reward():
                    return [
                        (
                            [(crit[x_j + y_j] - crit[x_j]) for crit in c_j]
                            if j in E else [0] * len(c_j)
                        )
                        for j, (c_j, y_j, x_j) in enumerate(zip(c, y, x))
                    ]
            else:
                # each job is an item with one copy (0-1 MDKP)
                n = len(jobs)
                x = [0] * n  # selected jobs
                E = set(range(n))  # jobs still free to select
                u = [1] * n
                a = list(map(self.job_weight, jobs))  # resource usage of jobs
                c = list(map(self.job_reward, jobs))  # job rewards

                def calc_reward():
                    return [c_j * y_j for c_j, y_j in zip(c, y)]

            b = [self.resources[name] for name in self.workflow.global_resources]  # resource capacities

            while True:
                # Step 2: compute effective capacities
                y = [
                        (
                            min(
                                (min(u[j], b_i // a_j_i) if a_j_i > 0 else u[j])
                                for b_i, a_j_i in zip(b, a[j]) if a_j_i
                            )
                            if j in E else 0
                        )
                        for j in range(n)
                    ]
                if not any(y):
                    break
                y = [(max(1, int(self.greedyness * y_j)) if y_j > 0 else 0) for y_j in y]

                # Step 3: compute rewards on cumulative sums
                reward = calc_reward()
                j_sel = max(E, key=reward.__getitem__)  # argmax

                # Step 4: batch increment
                y_sel = y[j_sel]

                # Step 5: update information
                x[j_sel] += y_sel
                b = [b_i - (a_j_i * y_sel) for b_i, a_j_i in zip(b, a[j_sel])]
                u[j_sel] -= y_sel
                if not u[j_sel] or self.greedyness == 1:
                    E.remove(j_sel)
                if not E:
                    break

            if self.select_by_rule:
                # Solution is the list of jobs that was selected from the selected rules
                solution = list(chain(
                    *[jobs[rules[j]][:x_] for j, x_ in enumerate(x)]))
            else:
                solution = [job for job, sel in zip(jobs, x) if sel]
            # update resources
            for name, b_i in zip(self.resources, b):
                self.resources[name] = b_i
            return solution

    def calc_resource(self, name, value):
        return min(value, self.workflow.global_resources[name])

    def rule_weight(self, rule):
        res = rule.resources
        return [
            self.calc_resource(name, res.get(name, 0))
            for name in self.workflow.global_resources]

    def rule_reward(self, rule, jobs=None):
        jobs = jobs[rule]
        return (
            self.priority_reward(jobs),
            self.downstream_reward(jobs),
            cumsum([job.inputsize for job in jobs]))

    def dryrun_rule_reward(self, rule, jobs=None):
        jobs = jobs[rule]
        return (
            self.priority_reward(jobs),
            self.downstream_reward(jobs),
            [0] * (len(jobs) + 1))

    def priority_reward(self, jobs):
        return cumsum(self.dag.priorities(jobs))

    def downstream_reward(self, jobs):
        return cumsum(self.dag.downstream_sizes(jobs))

    def thread_reward(self, jobs):
        """ Thread-based reward for jobs. Using this maximizes core 
        saturation, but does not lead to faster computation in general."""
        return cumsum([job.threads for job in jobs])

    def job_weight(self, job):
        res = job.resources_dict
        return [
            self.calc_resource(name, res.get(name, 0))
            for name in self.workflow.global_resources]

    def job_reward(self, job):
        return (
            self.dag.priority(job),
            self.dag.downstream_size(job),
            job.inputsize
        )

    def dryrun_job_reward(self, job):
        return (
            self.dag.priority(job),
            self.dag.downstream_size(job)
        )

    def progress(self):
        """ Display the progress. """
        logger.progress(done=self.finished_jobs, total=len(self.dag))
