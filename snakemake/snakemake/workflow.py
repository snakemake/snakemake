# -*- coding: utf-8 -*-

import re, sys, os, traceback, glob, signal
from multiprocessing import Event
from collections import defaultdict, OrderedDict
from tempfile import TemporaryFile

from snakemake.rules import Rule
from snakemake.exceptions import MissingOutputException, MissingInputException, AmbiguousRuleException, CyclicGraphException, MissingRuleException, RuleException, CreateRuleException, ProtectedOutputException, UnknownRuleException, NoRulesException
from snakemake.shell import shell, format
from snakemake.jobs import Job, KnapsackJobScheduler, ClusterJobScheduler
from snakemake.parser import compile_to_python
from snakemake.io import protected, temp

__author__ = "Johannes Köster"

class Jobcounter:
	def __init__(self):
		self._count = 0
		self._done = 0
	
	def add(self):
		self._count += 1
	
	def done(self):
		self._done += 1
	
	def __str__(self):
		return "{} of {} steps ({}%) done".format(self._done, self._count, int(self._done / self._count * 100))

class JobCounterSemaphore:
	def __init__(self, value):
		self.value = value
		self.event = Event()
	
	def release(self):
		self.value -= 1
		if self.value == 0:
			self.event.set()
	
	def wait(self):
		self.event.wait()

class Workflow:
	def __init__(self):
		"""
		Create the controller.
		"""
		self.init()
	
	def init(self, clear = False):
		if clear:
			for k in list(globals().keys()):
				if k not in self._virgin_globals:
					del globals()[k]
		else:
			self._virgin_globals = None
		self.__rules = OrderedDict()
		self.__last = None
		self.__first = None
		self.__altfirst = None
		self.__workdir_set = False
		self._jobs_finished = None
		self._runtimes = defaultdict(list)
		self._cores = 1
		self.rowmaps = dict()
		self.jobcounter = None
		self.rule_count = 0
		self.errors = False
	
	def get_cores(self):
		return self._cores
	
	def set_cores(self, cores):
		self._cores = cores
	
	def report_runtime(self, rule, runtime):
		self._runtimes[rule].append(runtime)
		
	def get_runtimes(self):
		for rule, runtimes in self._runtimes.items():
			s = sum(runtimes)
			yield rule, min(runtimes), max(runtimes), s, s / len(runtimes)

	def clear(self):
		self.init(clear = True)
		
	def set_job_finished(self, job = None, error = False):
		if error:
			self.errors = True
		self._jobs_finished.release()
		
	def get_snakefile_globals(self):
		return self.__snakefile_globals
	
	def get_pool(self):
		"""
		Return the current thread pool.
		"""
		return self.__pool
	
	def add_rule(self, name, lineno = None, snakefile = None):
		"""
		Add a rule.
		"""
		if self.is_rule(name):
			raise CreateRuleException("The name {} is already used by another rule".format(name))
		if "__" + name in globals():
			raise CreateRuleException("The name __{} is already used by a variable.".format(name))
		rule = Rule(name, self, lineno = lineno, snakefile = snakefile)
		self.__rules[rule.name] = rule
		self.__last = rule
		if not self.__first:
			self.__first = rule.name
			
	def is_rule(self, name):
		"""
		Return True if name is the name of a rule.
		
		Arguments
		name -- a name
		"""
		return name in self.__rules
	
	def has_run(self, rule):
		return "__" + rule.name in  globals()
	
	def get_run(self, rule):
		return globals()["__" + rule.name]
	
	def get_producers(self, files, exclude = None):
		for rule in self.get_rules():
			if rule != exclude:
				for f in files:
					if rule.is_producer(f):
						yield rule, f

	def get_rule(self, name):
		"""
		Get rule by name.
		
		Arguments
		name -- the name of the rule
		"""
		if not self.__rules:
			raise NoRulesException()
		if not name in self.__rules:
			raise UnknownRuleException(name)
		return self.__rules[name]

	def last_rule(self):
		"""
		Return the last rule.
		"""
		return self.__last

	def run_first_rule(self, dryrun = False, touch = False, forcethis = False, forceall = False, give_reason = False, cluster = None):
		"""
		Apply the rule defined first.
		"""
		first = self.__first
		if not first:
			for key, value in self.__rules.items():
				first = key
				break
		return self._run([(self.get_rule(first), None)], dryrun = dryrun, touch = touch, forcethis = forcethis, forceall = forceall, give_reason = give_reason, cluster = cluster)
			
	def get_file_producers(self, files, dryrun = False, forcethis = False, forceall = False):
		"""
		Return a dict of rules with requested files such that the requested files are produced.
		
		Arguments
		files -- the paths of the files to produce
		"""
		producers = dict()
		missing_input_ex = defaultdict(list)
		for rule, file in self.get_producers(files):
			try:
				rule.run(file, jobs=dict(), forceall = forceall, dryrun = True, visited = set())
				if file in producers:
					raise AmbiguousRuleException(producers[file], rule)
				producers[file] = rule
			except MissingInputException as ex:
				missing_input_ex[file].append(ex)
		
		toraise = []
		for file in files:
			if not file in producers:
				if file in missing_input_ex:
					toraise += missing_input_ex[file]
				else:
					toraise.append(MissingRuleException(file))
		if toraise:
			raise RuleException(include = toraise)

		return [(rule, file) for file, rule in producers.items()]
	
	def run_rules(self, targets, dryrun = False, touch = False, forcethis = False, forceall = False, give_reason = False, cluster = None):
		ruletargets, filetargets = [], []
		for target in targets:
			if workflow.is_rule(target):
				ruletargets.append(target)
			else:
				filetargets.append(target)
		
		torun = self.get_file_producers(filetargets, forcethis = forcethis, forceall = forceall, dryrun = dryrun) + \
			[(self.get_rule(name), None) for name in ruletargets]
				
		return self._run(torun, dryrun = dryrun, touch = touch, forcethis = forcethis, forceall = forceall, give_reason = give_reason, cluster = cluster)
	
	def _run(self, torun, dryrun = False, touch = False, forcethis = False, forceall = False, give_reason = False, cluster = None):
		self.jobcounter = Jobcounter()
		jobs = dict()
		
		for rule, requested_output in torun:
			job = rule.run(requested_output, jobs=jobs, forcethis = forcethis, forceall = forceall, dryrun = dryrun, give_reason = give_reason, touch = touch, visited = set(), jobcounter = self.jobcounter)
			job.add_callback(self.set_job_finished)

		self._jobs_finished = JobCounterSemaphore(len(torun))
		
		
		if cluster:
			scheduler = ClusterJobScheduler(set(jobs.values()), self, submitcmd = cluster)
		else:
			scheduler = KnapsackJobScheduler(set(jobs.values()), self)
		scheduler.schedule()


		self._jobs_finished.wait()
		scheduler.terminate()
		if self.errors:
			return 1
		return 0

	def check_rules(self):
		"""
		Check all rules.
		"""
		for rule in self.get_rules():
			rule.check()

	def get_rules(self):
		"""
		Get the list of rules.
		"""
		return self.__rules.values()

	def is_produced(self, files):
		"""
		Return True if files are already produced.
		
		Arguments
		files -- files to check
		"""
		for f in files:
			if not os.path.exists(f): return False
		return True
	
	def is_newer(self, files, time):
		"""
		Return True if files are newer than a time
		
		Arguments
		files -- files to check
		time -- a time
		"""
		for f in files:
			if os.stat(f).st_mtime > time: return True
		return False

	def include(self, snakefile, overwrite_first_rule = False):
		"""
		Include a snakefile.
		"""
		first_rule = self.__first
		code, rowmap, rule_count = compile_to_python(snakefile, rule_count = self.rule_count)
		self.rule_count += rule_count
		self.rowmaps[snakefile] = rowmap
		exec(compile(code, snakefile, "exec"), globals())
		if not overwrite_first_rule:
			self.__first = first_rule

	def set_workdir(self, workdir):
		if not self.__workdir_set:
			if not os.path.exists(workdir):
				os.makedirs(workdir)
			os.chdir(workdir)
			self.__workdir_set = True

workflow = Workflow()

def _include(path):
	workflow.include(path)

def _set_workdir(path):
	workflow.set_workdir(path)

def _add_rule(name, lineno = None, snakefile = None):
	workflow.add_rule(name, lineno = lineno, snakefile = snakefile)

def _set_input(*paths, **kwpaths):
	workflow.last_rule().set_input(*paths, **kwpaths)

def _set_output(*paths, **kwpaths):
	workflow.last_rule().set_output(*paths, **kwpaths)

def _set_message(message):
	workflow.last_rule().set_message(message)
	
def _set_threads(threads):
	workflow.last_rule().set_threads(threads)

workflow._virgin_globals = dict(globals())
