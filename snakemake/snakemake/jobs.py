# -*- coding: utf-8 -*-

import signal
import sys, time, os, threading, multiprocessing
from itertools import chain
from snakemake.exceptions import TerminatedException, MissingOutputException, RuleException, \
	ClusterJobException, print_exception
from snakemake.shell import shell
from snakemake.io import IOFile, temp, protected
from snakemake.logging import logger
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

if os.name == "posix":
	from multiprocessing import Event
	PoolExecutor = ProcessPoolExecutor
else:
	from threading import Event
	PoolExecutor = ThreadPoolExecutor

__author__ = "Johannes Köster"

def run_wrapper(run, rulename, ruledesc, input, output, wildcards, 
		threads, rowmaps, rulelineno, rulesnakefile):
	"""
	Wrapper around the run method that handles directory creation and
	output file deletion on error.
	
	Arguments
	run -- the run method
	input -- list of input files
	output -- list of output files
	wildcards -- so far processed wildcards
	"""
	logger.info(ruledesc)

	for o in output:
		o.prepare()
	try:
		t0 = time.time()
		# execute the actual run method.
		run(input, output, wildcards, threads)
		# finish all spawned shells.
		shell.join_all()
		runtime = time.time() - t0
		for o in output:
			o.created(rulename, rulelineno, rulesnakefile)
		for i in input:
			i.used()
		return runtime
	except (Exception, BaseException) as ex:
		# Remove produced output on exception
		for o in output:
			o.remove()
		if not isinstance(ex, TerminatedException):
			print_exception(ex, rowmaps)
			raise Exception()

class Job:
	count = 0

	@staticmethod
	def cleanup_unfinished(jobs):
		for job in jobs:
			job.cleanup()

	def __init__(self, workflow, rule = None, message = None, 
			input = None, output = None, wildcards = None, 
			threads = 1, depends = set(), dryrun = False, 
			touch = False, needrun = True):
		self.workflow = workflow
		self.rule = rule
		self.message = message
		self.input = input
		self.output = output
		self.wildcards = wildcards
		self.threads = threads
		self.dryrun = dryrun
		self.touch = touch
		self.needrun = needrun
		self.depends = set(depends)
		self.depending = list()
		self.is_finished = False
		self._callbacks = list()
		self._error_callbacks = list()
		self.jobid = Job.count
		Job.count += 1
		for other in self.depends:
			other.depending.append(self)

	def all_jobs(self):
		yield self
		for job in self.depends:
			for j in job.all_jobs():
				yield j
	
	def print_message(self):
		logger.info(self.message)
		
	def run(self, run_func):
		if not self.needrun:
			self.finished()
		elif self.dryrun:
			self.print_message()
			self.finished()
		elif self.touch:
			logger.info(self.message)
			for o in self.output:
				o.touch(self.rule.name, self.rule.lineno, self.rule.snakefile)
			# sleep shortly to ensure that output files of different rules 
			# are not touched at the same time.
			time.sleep(0.1)
			self.finished()
		else:
			run_func(self)
	
	def get_run_args(self):
		return (self.rule.get_run(), self.rule.name, self.message, 
			self.input, self.output, self.wildcards, self.threads, 
			self.workflow.rowmaps, self.rule.lineno, self.rule.snakefile)
	
	def add_callback(self, callback):
		""" Add a callback that is invoked when job is finished. """
		self._callbacks.append(callback)

	def add_error_callback(self, callback):
		self._error_callbacks.append(callback)
	
	def finished(self, future = None):
		""" Set job to be finished. """	
		if future:
			ex = future.exception()
			if ex:
				for callback in self._error_callbacks:
					callback()
				return
		self.is_finished = True
		if self.needrun and not self.dryrun:
			self.workflow.jobcounter.done()
			logger.info(self.workflow.jobcounter)
			if future != None:
				self.workflow.report_runtime(self.rule, future.result())
		for other in self.depending:
			other.depends.remove(self)
		for callback in self._callbacks:
			callback(self)

	def cleanup(self):
		if not self.is_finished:
			for o in self.output:
				o.remove()

	def dot(self):
		label = self.rule.name
		new_wildcards = self.new_wildcards()
		if not self.depends or new_wildcards:
			for wildcard, value in new_wildcards:
				label += "\\n{}: {}".format(wildcard, value)
		edges = ("{} -> {};".format(j.jobid, self.jobid) 
			for j in self.depends if j.needrun)
		node = ('{}[label = "{}"];'.format(self.jobid, label),)
		return chain(node, edges)

	def new_wildcards(self):
		new_wildcards = set(self.wildcards.items())
		for job in self.depends:
			if not new_wildcards:
				return set()
			for wildcard in job.wildcards.items():
				new_wildcards.discard(wildcard)
		return new_wildcards

	def __repr__(self):
		return self.rule.name


class KnapsackJobScheduler:
	def __init__(self, jobs, workflow):
		""" Create a new instance of KnapsackJobScheduler. """
		self.workflow = workflow
		self._maxcores = workflow.cores if workflow.cores else multiprocessing.cpu_count()
		self._cores = self._maxcores
		self._pool = PoolExecutor(max_workers = self._cores)
		self._jobs = set(jobs)
		self._open_jobs = Event()
		self._open_jobs.set()
		self._errors = False

	def schedule(self):
		""" Schedule jobs that are ready, maximizing cpu usage. """
		while True:
			self._open_jobs.wait()
			self._open_jobs.clear()
			if self._errors:
				logger.warning("Will exit after finishing currently running jobs.")
				self._pool.shutdown()
				return False
			if not self._jobs:
				self._pool.shutdown()
				return True

			needrun, norun = [], set()
			for job in self._jobs:
				if job.depends:
					continue
				if job.needrun:
					if job.threads > self._maxcores:
						# reduce the number of threads so that it 
						# fits to available cores.
						if not job.dryrun:
							logger.warn(
								"Rule {} defines too many threads ({}), Scaling down to {}."
								.format(job.rule, job.threads, self._maxcores))
						job.threads = self._maxcores
					needrun.append(job)
				else: norun.add(job)
			
			run = self._knapsack(needrun)
			self._jobs -= run
			self._jobs -= norun
			self._cores -= sum(job.threads for job in run)
			for job in chain(run, norun):
				job.add_callback(self._finished)
				job.add_error_callback(self._error)
				job.run(self._run_job)
			
	
	def _run_job(self, job):
		future = self._pool.submit(run_wrapper, *job.get_run_args())
		future.add_done_callback(job.finished)
		
	def _finished(self, job):
		if job.needrun:
			self._cores += job.threads
		self._open_jobs.set()
	
	def _error(self):
		# clear jobs and stop the workflow
		self._errors = True
		self._jobs = set()
		self._open_jobs.set()
	
	def _knapsack(self, jobs):
		""" Solve 0-1 knapsack to maximize cpu utilization. """
		dimi, dimj = len(jobs) + 1, self._cores + 1
		K = [[0 for c in range(dimj)] for i in range(dimi)]
		for i in range(1, dimi):
			for j in range(1, dimj):
				t = jobs[i-1].threads
				if t > j:
					K[i][j] = K[i - 1][j]
				else:
					K[i][j] = max(K[i - 1][j], t + K[i - 1][j - t])
		
		solution = set()
		i = dimi - 1
		j = dimj - 1
		while i > 0:
			if K[i][j] != K[i-1][j]:
				job = jobs[i - 1]
				solution.add(job)
				j = j - job.threads
			i -= 1
		
		return solution

class ClusterJobScheduler:
	def __init__(self, jobs, workflow, submitcmd = "qsub"):
		self.workflow = workflow
		self._jobs = set(jobs)
		self._submitcmd = submitcmd
		self._open_jobs = Event()
		self._open_jobs.set()
		self._error = False
		self._cores = workflow.cores

	def schedule(self):
		while True:
			self._open_jobs.wait()
			self._open_jobs.clear()
			if self._error:
				logger.warning("Will exit after finishing currently running jobs.")
				return False
			if not self._jobs:
				return True
			needrun, norun = set(), set()
			for job in self._jobs:
				if job.depends:
					continue
				if job.needrun:
					needrun.add(job)
				else: norun.add(job)

			self._jobs -= needrun
			self._jobs -= norun
			for job in chain(needrun, norun):
				job.add_callback(self._finished)
				job.run(self._run_job)
	
	def _run_job(self, job):
		job.print_message()
		workdir = os.getcwd()
		prefix = ".snakemake"
		jobid = "_".join(job.output).replace("/", "_")
		jobscript = "{}.{}.sh".format(prefix, jobid)
		jobfinished = "{}.{}.jobfinished".format(prefix, jobid)
		jobfailed = "{}.{}.jobfailed".format(prefix, jobid)
		cores = self._cores if self._cores else ""
		shell("""
			echo '#!/bin/sh' > "{jobscript}"
			echo 'snakemake --force -j{self._cores} --directory {workdir} --nocolor --quiet {job.output} && touch "{jobfinished}" || touch "{jobfailed}"' >> "{jobscript}"
			chmod +x "{jobscript}"
			{self._submitcmd} "{jobscript}"
		""")
		threading.Thread(target=self._wait_for_job, args=(job, jobfinished, jobfailed, jobscript)).start()

	def _finished(self, job):
		self._open_jobs.set()
		
	def _wait_for_job(self, job, jobfinished, jobfailed, jobscript):
		while True:
			if os.path.exists(jobfinished):
				os.remove(jobfinished)
				os.remove(jobscript)
				job.finished()
				return
			if os.path.exists(jobfailed):
				os.remove(jobfailed)
				os.remove(jobscript)
				print_exception(ClusterJobException(job), self.workflow.rowmaps)
				self._error = True
				self._jobs = set()
				self._open_jobs.set()
				return
			time.sleep(1)

def print_job_dag(jobs):
	print("digraph snakemake_dag {")
	for job in jobs:
		for edge in job.dot():
			print("\t" + edge)
	print("}")
