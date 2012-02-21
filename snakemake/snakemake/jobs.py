import sys, os, time, stat, traceback
from snakemake.exceptions import MissingOutputException, RuleException, print_exception
from snakemake.shell import shell

class IOFile(str):
	_register = dict()
	
	@classmethod
	def create(cls, file, temp = False):
		obj = None
		if file in cls._register:
			obj = cls._register[file]
		else:
			obj = IOFile(file)
			print(obj)
			cls._register[file] = obj 
		if obj._temp == None:
			obj._temp = temp
		return obj
	
	def __init__(self, file, protected = False):
		self._file = file
		self._needed = 0
		self._temp = None
		self._protected = protected
		
	def need(self):
		self._needed += 1
	
	def used(self):
		self._needed -= 1
		if self._temp:
			os.remove(self._file)
	
	def prepare(self):
		dir = os.path.dirname(self._file)
		if len(dir) > 0 and not os.path.exists(dir):
			os.makedirs(dir)
	
	def created(self, rulename, lineno, snakefile):
		if not os.path.exists(self._file):
			raise MissingOutputException("Output file {} not produced by rule {}.".format(self._file, rulename), lineno = lineno, snakefile = snakefile)
		if self._protected:
			mode = os.stat(o).st_mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH
			if os.path.isdir(o):
				for root, dirs, files in os.walk(o):
					for d in dirs:
						os.chmod(os.path.join(o, d), mode)
					for f in files:
						os.chmod(os.path.join(o, f), mode)
			else:
				os.chmod(o, mode)
	
	def remove(self):
		if os.path.exists(self._file):
			if os.path.isdir(self._file):
				for root, dirs, files in os.walk(o):
					for f in files:
						os.remove(f)
				for root, dirs, files in os.walk(o):
					for d in dirs:
						os.rmdir(d)
			else:
				os.remove(o)

	def apply_wildcards(self, wildcards):
		return self.create(self._file.format(**wildcards))
	
	def __str__(self):
		return self._file

def temp(file):
	return IOFile.create(file, temp = True)

def protected(file):
	return IOFile.create(file, protected = True)

def run_wrapper(run, rulename, ruledesc, input, output, wildcards, rowmaps, rulelineno, rulesnakefile):
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
		run(input, output, wildcards)
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
	def __init__(self, workflow, rule = None, message = None, input = None, output = None, wildcards = None, depends = set(), dryrun = False, needrun = True):
		self.rule, self.message, self.input, self.output, self.wildcards = rule, message, input, output, wildcards
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
				self.workflow.get_pool().apply_async(
					run_wrapper, 
					(self.rule.get_run(), self.rule.name, self.message, self.input, self.output, self.wildcards, self.workflow.rowmaps, self.rule.lineno, self.rule.snakefile), 
					callback=self._finished_callback, 
					error_callback=self._raise_error
				)
		else:
			self._wakeup_waiting()
	
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
