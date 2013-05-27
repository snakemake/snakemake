# -*- coding: utf-8 -*-

import os
import sys
import threading
import multiprocessing
import operator
import signal
from functools import partial

from snakemake.executors import DryrunExecutor, TouchExecutor
from snakemake.executors import ClusterExecutor, CPUExecutor
from snakemake.logging import logger

__author__ = "Johannes Köster"


class JobScheduler:
    def __init__(
        self,
        workflow,
        dag,
        cores,
        dryrun=False,
        touch=False,
        cluster=None,
        immediate_submit=False,
        quiet=False,
        printreason=False,
        printshellcmds=False,
        keepgoing=False,
        output_wait=3):
        """ Create a new instance of KnapsackJobScheduler. """
        self.dag = dag
        self.workflow = workflow
        self.dryrun = dryrun
        self.quiet = quiet
        self.keepgoing = keepgoing
        self.maxcores = cores
        self.running = set()
        self.finished_jobs = 0
        self._cores = self.maxcores
        use_threads = os.name != "posix"
        self._open_jobs = (multiprocessing.Event() if not use_threads
            else threading.Event())
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
                output_wait=output_wait)
        elif touch:
            self._executor = TouchExecutor(
                workflow, dag, printreason=printreason,
                quiet=quiet, printshellcmds=printshellcmds,
                output_wait=output_wait)
        elif cluster:
            # TODO properly set cores
            self._executor = ClusterExecutor(
                workflow, dag, None, submitcmd=cluster,
                printreason=printreason, quiet=quiet,
                printshellcmds=printshellcmds, output_wait=output_wait)
            self._open_jobs = threading.Event()
            self._job_weight = self.simple_job_weight
            if immediate_submit:
                self._submit_callback = partial(
                    self._proceed,
                    update_dynamic=False,
                    print_progress=False,
                    update_resources=False)
        else:
            self._executor = CPUExecutor(
                workflow, dag, cores, printreason=printreason,
                quiet=quiet, printshellcmds=printshellcmds,
                threads=use_threads,
                output_wait=output_wait)
        self._open_jobs.set()

    @property
    def stats(self):
        try:
            return self._executor.stats
        except AttributeError:
            raise TypeError("Executor does not support stats")

    def candidate(self, job):
        """ Return whether a job is a candidate to be executed. """
        return (job not in self.running and (self.dryrun or 
            (not job.dynamic_input and not self.dag.dynamic(job))))

    @property
    def open_jobs(self):
        """ Return open jobs. """
        return filter(self.candidate, list(self.dag.ready_jobs))

    def schedule(self):
        """ Schedule jobs that are ready, maximizing cpu usage. """
        while True:
            try:
                self._open_jobs.wait()
            except:
                # this will be caused becaus of SIGTERM or SIGINT, so exit with error
                self._errors = True
            self._open_jobs.clear()
            if not self.keepgoing and self._errors:
                logger.warning("Will exit after finishing "
                    "currently running jobs.")
                self._executor.shutdown()
                return False
            if not any(self.open_jobs):
                self._executor.shutdown()
                return not self._errors

            needrun = list()
            for job in self.open_jobs:
                if job.threads > self.maxcores:
                    job.threads = self.maxcores
                needrun.append(job)
            assert needrun

            logger.debug("Ready jobs:\n\t" + "\n\t".join(map(str, needrun)))

            run = self.job_selector(needrun)
            logger.debug("Selected jobs:\n\t" + "\n\t".join(map(str, run)))
            self.running.update(run)
            self._cores -= sum(job.threads for job in run)
            for job in run:
                self._executor.run(
                    job, callback=self._finish_callback,
                    submit_callback=self._submit_callback,
                    error_callback=self._error)

    def _noop(self, job):
        pass

    def _proceed(
        self, job, update_dynamic=True, print_progress=False,
        update_resources=True):
        """ Do stuff after job is finished. """
        if update_resources:
            self.finished_jobs += 1
            self.running.remove(job)
            self._cores += job.threads

        self.dag.finish(job, update_dynamic=update_dynamic)

        if print_progress:
            self.progress()

        if any(self.open_jobs) or not self.running:
            # go on scheduling if open jobs are ready or no job is running any
            self._open_jobs.set()

    def _error(self):
        """ Clear jobs and stop the workflow. """
        self._errors = True
        if self.keepgoing:
            logger.warning("Job failed, going on with independent jobs.")
        else:
            self._open_jobs.set()

    def job_selector(self, jobs):
        """ Solve 0-1 knapsack to maximize cpu utilization. """
        resources = self.workflow.resources

        dimi, dimj = len(jobs) + 1, self._cores + 1
        K = [[(0, 0) for c in range(dimj)] for i in range(dimi)]
        for i in range(1, dimi):
            for j in range(1, dimj):
                job = jobs[i - 1]
                w = self.job_weight(job)
                v = job.priority, job.inputsize if not self.dryrun else 0
                if w > j:
                    K[i][j] = K[i - 1][j]
                else:
                    K[i][j] = max(K[i - 1][j],
                        tuple(map(operator.add, v, K[i - 1][j - w])))

        solution = set()
        i = dimi - 1
        j = dimj - 1
        while i > 0:
            if K[i][j] != K[i - 1][j]:
                job = jobs[i - 1]
                solution.add(job)
                j = j - job.threads
            i -= 1
        return solution

    def job_selector_greedy(self, jobs):
        """
        Using the greedy heuristic from 
        "A Greedy Algorithm for the General Multidimensional Knapsack
Problem", Akcay, Li, Xu, Annals of Operations Research, 2012
        """
        # solve over the rules instead of jobs (much less)
        _jobs = defaultdict(list)
        _jobs.update((job.rule, job) for job in jobs)
        jobs = _jobs
        rules = list(jobs)

        # Step 1: initialization
        x = [0] * len(rules)
        E = set(range(len(rules)))
        u = [len(jobs[rule]) for rule in rules]
        b = [self._cores] + [val for val in self.workflow.resources.values()]
        a = list(map(self.rule_weight, rules))
        c = list(map(self.rule_reward, rules))

        while True:
            # Step 2: compute effective capacities
            y = [
                min(b_i // a_j_i for b_i, a_j_i in zip(b, a[j]) if a_j_i)
                for j in E]
            if not any(y):
                break

            # Step 3: compute rewards
            reward = [c_j * y_j for c_j, y_j in zip(c, y)]
            j_ = max(E, key=reward.__getitem__) # argmax

            # Step 4: batch increment
            y_ = min(u[j_], max(1, alpha * y[j_]))
            
            # Step 5: update information
            x[j_] += y_
            b = [b_i - (a_j_i * y_) for b_i, a_j_i in zip(b, a[j_])]
            u[j_] -= y_
            if not u[j_] or alpha == 1:
                E.remove(j_)
            if not E:
                break
                
        # Solution is the list of jobs that was selected from the selected rules
        solution = list(chain(*[jobs[rules[j]][:x_] for j, x_ in enumerate(x)]))
        return solution

    def job_weight(self, job):
        """ Job weight that uses threads. """
        return job.threads

    def simple_job_weight(self, job):
        """ Job weight that ignores threads. """
        return 1

    def progress(self):
        """ Display the progress. """
        logger.info("{} of {} steps ({:.0%}) done".format(self.finished_jobs,
            len(self.dag), self.finished_jobs / len(self.dag)))
