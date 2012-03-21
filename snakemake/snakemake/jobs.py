import sys, time
from threading import Thread
from snakemake.exceptions import MissingOutputException, RuleException, print_exception
from snakemake.shell import shell
from snakemake.io import IOFile, temp, protected
from multiprocessing import Process, Pool
from itertools import chain

def run_wrapper(run, rulename, ruledesc, input, output, wildcards, threads, rowmaps, rulelineno, rulesnakefile):
	"""
	Wrapper around the run method that handles directory creation and output file deletion on error.
	
	Arguments
	run -- the run method
	input -- list of input files
	output -- list of output files
	wildcards -- so far processed wildcards
	"""
	print(ruledesc)

	for o in output:
		o.prepare()
		
	try:
		t0 = time.time()
		# execute the actual run method.
		run(input, output, wildcards, threads)
		# finish all spawned shells.
		shell.join()
		runtime = time.time() - t0
		for o in output:
			o.created(rulename, rulelineno, rulesnakefile)
		for i in input:
			i.used()
		return runtime
	except (Exception, BaseException) as ex:
		print_exception(ex, rowmaps)
		
		# Remove produced output on exception
		for o in output:
			o.remove()
		raise Exception()

class Job:
	def __init__(self, workflow, rule = None, message = None, input = None, output = None, wildcards = None, threads = 1, depends = set(), dryrun = False, needrun = True):
		self.workflow = workflow
		self.rule = rule
		self.message = message
		self.input = input
		self.output = output
		self.wildcards = wildcards
		self.threads = threads
		self.dryrun = dryrun
		self.needrun = needrun
		self.depends = depends
		self.depending = list()
		self._callbacks = list()
		for other in self.depends:
			other.depending.append(self)
	
	def add_callback(self, callback):
		""" Add a callback that is invoked when job is finished. """
		self._callbacks.append(callback)
	
	def finished(self, runtime = None):
		""" Set job to be finished. """
		if self.needrun:
			self.workflow.jobcounter.done()
			print(self.workflow.jobcounter)
			if runtime != None:
				self.workflow.report_runtime(self.rule, runtime)
		for callback in self._callbacks:
			callback(self)


class KnapsackJobScheduler:
	def __init__(self, jobs, workflow):
		""" Create a new instance of KnapsackJobScheduler. """
		self.workflow = workflow
		self._cores = workflow.get_cores()
		self._pool = Pool(self._cores)
		self._jobs = set(jobs)
	
	def schedule(self):
		""" Schedule jobs that are ready, maximizing cpu usage. """
		needrun, norun = [], set()
		for job in self._jobs:
			if not job.depends:
				if job.threads > self._cores:
					# reduce the number of threads so that it fits to available cores.
					job.threads = self._cores
				if job.needrun: needrun.append(job)
				else: norun.add(job)
		
		run = self._knapsack(needrun)
		self._jobs -= run
		self._jobs -= norun
		self._cores -= sum(job.threads for job in run)
		for job in chain(run, norun):
			def finished(runtime = None):
				""" Schedule remaining jobs. """
				job.finished(runtime)
				# remove this job from others dependencies
				for other in job.depending:
					other.depends.remove(job)
				self._cores += job.threads
				self.schedule()
			
			if not job.needrun:
				finished()
			elif not job.dryrun:
				self._pool.apply_async(
					run_wrapper, 
					(job.rule.get_run(), job.rule.name, job.message, job.input, job.output, job.wildcards, job.threads, job.workflow.rowmaps, job.rule.lineno, job.rule.snakefile),
					callback = finished,
					error_callback = self._error
				)
			else:
				print(job.message)
				finished()
	
	def _error(self, error):
		# simply stop because exception was printed in run_wrapper
		self.workflow.set_job_finished(error = True)
	
	def _knapsack(self, jobs):
		""" Solve 0-1 knapsack to maximize cpu utilization. """
		K = [[0 for c in range(self._cores + 1)] for i in range(len(jobs) + 1)]
		for i in range(1, len(jobs) + 1):
			for j in range(1, self._cores + 1):
				t = jobs[i-1].threads
				if t > j:
					K[i][j] = K[i - 1][j]
				else:
					K[i][j] = max(K[i - 1][j], t + K[i - 1][j - t])
		
		solution = set()
		i = len(jobs)
		j = self._cores
		while i > 0:
			if K[i][j] != K[i-1][j]:
				job = jobs[i - 1]
				solution.add(job)
				j = j - job.threads
			i -= 1
		
		return solution

