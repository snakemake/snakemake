
import os, threading
from itertools import filterfalse

from snakemake.executors import DryrunExecutor, TouchExecutor, ClusterExecutor, CPUExecutor
from snakemake.stats import Stats
from snakemake.logging import logger

class JobScheduler:
	def __init__(self, workflow, dag, cores, dryrun = False, touch = False, cluster = None, quiet = False, printreason = False, printshellcmds = False, keepgoing = False):
		""" Create a new instance of KnapsackJobScheduler. """
		self.dag = dag
		self.dryrun = dryrun
		self.quiet = quiet
		self.keepgoing = keepgoing
		self.maxcores = cores
		self.running = set()
		self.finished_jobs = 0
		self.stats = Stats()
		self._cores = self.maxcores
		use_threads = os.name == "posix"
		self._open_jobs = multiprocessing.Event() if not use_threads else threading.Event()
		self._errors = False
		self._finished = False
		self._job_queue = None
		
		if dryrun:
			self._executor = DryrunExecutor(workflow, dag, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
		elif touch:
			self._executor = TouchExecutor(workflow, dag, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
		elif cluster:
			# TODO properly set cores
			self._executor = ClusterExecutor(workflow, dag, None, submitcmd=cluster, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
			self._open_jobs = threading.Event()
			self._job_weight = self.simple_job_weight
		else:
			self._executor = CPUExecutor(workflow, dag, cores, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds, threads=use_threads)
		self._open_jobs.set()
	
	def candidate(self, job):
		return job not in self.running and not self.dag.dynamic(job) and (self.dryrun or not job.dynamic_input)
	
	@property
	def open_jobs(self):
		return filter(self.candidate, self.dag.ready_jobs)
	
	def schedule(self):
		""" Schedule jobs that are ready, maximizing cpu usage. """
		while True:
			self._open_jobs.wait()
			self._open_jobs.clear()
			if not self.keepgoing and self._errors:
				logger.warning("Will exit after finishing currently running jobs.")
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

			run = self.job_selector(needrun)
			self.running.update(run)
			self._cores -= sum(job.threads for job in run)
			for job in run:
				self.stats.report_job_start(job)
				self._executor.run(job, callback=self._finish_job, error_callback=self._error)
		
	def _finish_job(self, job):
		self.stats.report_job_end(job)
		self.finished_jobs += 1
		self.running.remove(job)
		self.dag.finish(job, update_dynamic=not self.dryrun)
		self._cores += job.threads
		#import pdb; pdb.set_trace()
		if not self.quiet and not self.dryrun:
			self.progress()
		if any(self.open_jobs) or not self.running:
			# go on scheduling if open jobs are ready or no job is running any more
			self._open_jobs.set()
	
	def _error(self):
		# clear jobs and stop the workflow
		self._errors = True
		if self.keepgoing:
			logger.warning("Job failed, going on with independent jobs.")
		else:
			self._open_jobs.set()
	
	def job_selector(self, jobs):
		""" Solve 0-1 knapsack to maximize cpu utilization. """
		dimi, dimj = len(jobs) + 1, self._cores + 1
		K = [[0 for c in range(dimj)] for i in range(dimi)]
		for i in range(1, dimi):
			for j in range(1, dimj):
				job = jobs[i-1]
				w = self.job_weight(job)
				v = job.priority
				if w > j:
					K[i][j] = K[i - 1][j]
				else:
					K[i][j] = max(K[i - 1][j], v + K[i - 1][j - w])
		
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
	
	def job_weight(self, job):
		return job.threads
	
	def simple_job_weight(self, job):
		return 1
	
	def progress(self):
		logger.info("{} of {} steps ({:.0%}) done".format(self.finished_jobs, len(self.dag), self.finished_jobs / len(self.dag)))
