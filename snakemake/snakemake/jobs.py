import signal
import sys, time, os
from threading import Thread
from snakemake.exceptions import TerminatedException, MissingOutputException, RuleException, print_exception
from snakemake.shell import shell
from snakemake.io import IOFile, temp, protected
from snakemake.logging import logger
from multiprocessing import Process, Pool, Lock
import threading
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
	def __init__(self, workflow, rule = None, message = None, input = None, output = None, wildcards = None, threads = 1, depends = set(), dryrun = False, touch = False, needrun = True):
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
		self._callbacks = list()
		for other in self.depends:
			other.depending.append(self)
			
	def run(self, run_func):
		if not self.needrun:
			self.finished()
		elif self.dryrun:
			logger.info(self.message)
			self.finished()
		elif self.touch:
			logger.info(self.message)
			for o in self.output:
				o.touch(self.rule.lineno, self.rule.snakefile)
			# sleep shortly to ensure that output files of different rules are not touched at the same time.
			time.sleep(0.1)
			self.finished()
		else:
			run_func(self)
	
	def get_run_args(self):
		return (self.rule.get_run(), self.rule.name, self.message, self.input, self.output, self.wildcards, self.threads, self.workflow.rowmaps, self.rule.lineno, self.rule.snakefile)
	
	def add_callback(self, callback):
		""" Add a callback that is invoked when job is finished. """
		self._callbacks.append(callback)
	
	def finished(self, runtime = None):
		""" Set job to be finished. """	
		if self.needrun and not self.dryrun:
			self.workflow.jobcounter.done()
			logger.info(self.workflow.jobcounter)
			if runtime != None:
				self.workflow.report_runtime(self.rule, runtime)
		for other in self.depending:
			other.depends.remove(self)
		for callback in self._callbacks:
			callback(self)
			
	def __repr__(self):
		return self.rule.name


class KnapsackJobScheduler:
	def __init__(self, jobs, workflow):
		""" Create a new instance of KnapsackJobScheduler. """
		self.workflow = workflow
		self._maxcores = workflow.get_cores()
		self._cores = self._maxcores
		self._pool = Pool(self._cores, maxtasksperchild = 1)
		self._jobs = set(jobs)
		self._lock = Lock()

	def terminate(self):
		self._pool.close()
		self._pool.terminate()
	
	def schedule(self):
		""" Schedule jobs that are ready, maximizing cpu usage. """
		self._lock.acquire()
		#import pdb; pdb.set_trace()
		needrun, norun = [], set()
		for job in self._jobs:
			if job.depends:
				continue
			if job.needrun:
				if job.threads > self._maxcores:
					# reduce the number of threads so that it fits to available cores.
					if not job.dryrun:
						logger.warn("Rule {} defines too many threads ({}), Scaling down to {}.".format(job.rule, job.threads, self._maxcores))
					job.threads = self._maxcores
				needrun.append(job)
			else: norun.add(job)
		
		run = self._knapsack(needrun)
		self._jobs -= run
		self._jobs -= norun
		self._cores -= sum(job.threads for job in run)
		self._lock.release()
		for job in chain(run, norun):
			job.add_callback(self._finished)
			job.run(self._run_job)
	
	def _run_job(self, job):
		self._pool.apply_async(
			run_wrapper, 
			job.get_run_args(),
			callback = job.finished,
			error_callback = self._error
		)
		
	def _finished(self, job):
		self._cores += job.threads
		self.schedule()
	
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

class ClusterJobScheduler:
	def __init__(self, jobs, workflow, submitcmd = "qsub"):
		self.workflow = workflow
		self._jobs = set(jobs)
		self._submitcmd = submitcmd
		self._lock = threading.Lock()

	def terminate(self):
		pass

	def schedule(self):
		self._lock.acquire()
		needrun, norun = set(), set()
		for job in self._jobs:
			if job.depends:
				continue
			print(job.rule)
			if job.needrun:
				needrun.add(job)
			else: norun.add(job)

		self._jobs -= needrun
		self._jobs -= norun
		self._lock.release()
		for job in chain(needrun, norun):
			job.add_callback(self._finished)
			job.run(self._run_job)
	
	def _run_job(self, job):
		prefix = ".snakemake"
		jobid = "{}{}".format("_".join(job.output), time.time())
		jobscript = "{}.{}.sh".format(prefix, jobid)
		jobguard = "{}.{}.jobguard".format(prefix, jobid)
		shell("""
			echo '#!/bin/sh' > {jobscript}
			echo 'snakemake {job.output} && touch {jobguard}' >> {jobscript}
			chmod +x {jobscript}
			{self._submitcmd} {jobscript}
		""")
		threading.Thread(target=self._wait_for_job, args=(job, jobguard, jobscript)).start()

	def _finished(self, job):
		self.schedule()
		
	def _wait_for_job(self, job, jobguard, jobscript):
		while not os.path.exists(jobguard):
			time.sleep(1)
		os.remove(jobguard)
		os.remove(jobscript)
		job.finished()
