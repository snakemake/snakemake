
import os, threading

from snakemake.executors import DryrunExecutor, TouchExecutor, ClusterExecutor, CPUExecutor
from snakemake.stats import Stats
from snakemake.logging import logger

class JobScheduler:
	def __init__(self, workflow, dag, cores, dryrun = False, touch = False, cluster = None, quiet = False, printreason = False, printshellcmds = False):
		""" Create a new instance of KnapsackJobScheduler. """
		self.dag = dag
		self.dryrun = dryrun
		self.quiet = quiet
		self.maxcores = cores
		self.finished_jobs = 0
		self.stats = Stats()
		self._cores = self.maxcores
		use_threads = os.name == "posix"
		self._open_jobs = multiprocessing.Event() if not use_threads else threading.Event()
		self._errors = False
		if dryrun:
			self._executor = DryrunExecutor(workflow, dag, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
			self.progress = lambda: None
		elif touch:
			self._executor = TouchExecutor(workflow, dag, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
		elif cluster:
			# TODO properly set cores
			self._executor = ClusterExecutor(workflow, dag, None, submitcmd=cluster, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
			self._open_jobs = threading.Event()
		else:
			self._executor = CPUExecutor(workflow, dag, cores, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds, threads=use_threads)
			self._selector = self._thread_based_selector
		self._open_jobs.set()

	def schedule(self):
		""" Schedule jobs that are ready, maximizing cpu usage. """
		while True:
			self._open_jobs.wait()
			self._open_jobs.clear()
			if self._errors:
				logger.warning("Will exit after finishing currently running jobs.")
				self._executor.shutdown()
				return False

			#import pdb; pdb.set_trace()
			needrun = list()
			for job in filter(self.dag.needrun, self.dag.ready_jobs):
				if job.threads > self.maxcores:
					# reduce the number of threads so that it 
					# fits to available cores.
					if not self.dryrun:
						logger.warn(
							"Rule {} defines too many threads ({}), Scaling down to {}."
							.format(job.rule, job.threads, self.maxcores))
					job.threads = self.maxcores
				needrun.append(job)
			if not needrun:
				self._executor.shutdown()
				return True

			run = self._selector(needrun)
			self._cores -= sum(job.threads for job in run)
			for job in run:
				self.stats.report_job_start(job)
				self._executor.run(job, callback=self._finished, error_callback=self._error)
		
	def _finished(self, job):
		self.stats.report_job_end(job)
		self.finished_jobs += 1
		needrun = self.dag.needrun(job)
		self.dag.finish(job, update_dynamic=not self.dryrun)
		self._cores += job.threads
		if not self.quiet:
			self.progress()
		self._open_jobs.set()
	
	def _error(self):
		# clear jobs and stop the workflow
		self._errors = True
		self._jobs = set()
		self._open_jobs.set()
	
	def _selector(self, jobs):
		return jobs[:self._cores]
	
	def _thread_based_selector(self, jobs):
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
	
	def progress(self):
		logger.info("{} of {} steps ({:.0%}) done".format(self.finished_jobs, len(self.dag), self.finished_jobs / len(self.dag)))
