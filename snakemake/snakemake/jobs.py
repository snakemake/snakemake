import sys, time
from threading import Thread
from snakemake.exceptions import MissingOutputException, RuleException, print_exception
from snakemake.shell import shell
from snakemake.io import IOFile, temp, protected
from multiprocessing import Process

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
		self.rule, self.message, self.input, self.output, self.wildcards, self.threads = rule, message, input, output, wildcards, threads
		self.dryrun, self.needrun = dryrun, needrun
		self.depends = depends
		self._callbacks = list()
		self.workflow = workflow
		self._queued = False
		self._finished = False
	
	def run(self, callback = None):
		if callback:
			if self._finished:
				callback(self)
				return
			self._callbacks.append(callback)
		if not self._queued:
			self._queued = True
			if self.depends:
				for job in list(self.depends):
					if job in self.depends:
						job.run(callback = self._wakeup_if_ready)
			else:
				self._wakeup()
	
	def _wakeup_if_ready(self, job):
		if job in self.depends:
			self.depends.remove(job)
		if not self.depends:
			self._wakeup()
	
	def _wakeup(self):
		if self.rule.has_run() and self.needrun:
			if self.dryrun:
				print(self.message)
				self._wakeup_waiting()
			else:
				self._apply()
		else:
			self._wakeup_waiting()
	
	def _apply(self):
		self.workflow.get_pool().apply_async(
			run_wrapper, 
			(self.rule.get_run(), self.rule.name, self.message, self.input, self.output, self.wildcards, self.threads, self.workflow.rowmaps, self.rule.lineno, self.rule.snakefile), 
			callback=self._finished_callback,
			error_callback=self._raise_error
		)
	
	def _finished_callback(self, value = None):
		self.workflow.jobcounter.done()
		print(self.workflow.jobcounter)
		if value != None:
			self.workflow.report_runtime(self.rule, value)
		self._wakeup_waiting()
		
	def _wakeup_waiting(self):
		self._finished = True
		for callback in self._callbacks:
			callback(self)
	
	def _raise_error(self, error):
		# simply stop because exception was printed in run_wrapper
		self.workflow.set_job_finished(error = True)

class JobScheduler:
	def __init__(self, jobs, cores, workflow):
		self.workflow = workflow
		self.jobs = jobs
		self.cores = cores
		self.pool = Pool()
	
	def schedule(self):
		ready = [job for job in self.jobs if not job.depends]
		run = self.knapsack(ready)
		for job in run:
			self.pool.apply_async(
				run_wrapper, 
				(job.rule.get_run(), job.rule.name, job.message, job.input, job.output, job.wildcards, job.threads, job.workflow.rowmaps, job.rule.lineno, job.rule.snakefile),
				callback = self.schedule,
				error_callback = self.error
			)

	
	def knapsack(self, jobs):
		R = [[0 for j in range(self.cores)] for i in range(len(jobs))]
		for i in range(len(jobs) - 1, -1, -1):
			for j in range(self.cores):
				if jobs[i].threads <= j + 1:
					R[i][j] = max(jobs[i].threads + R[i + 1][j - jobs[i].threads], R[i + 1][j])
		# backtrack best solution here
		return solution

	def error(self, error):
		# simply stop because exception was printed in run_wrapper
		self.workflow.set_job_finished(error = True)
		
