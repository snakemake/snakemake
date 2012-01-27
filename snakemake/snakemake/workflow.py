# -*- coding: utf-8 -*-

'''
Created on 13.11.2011

@author: Johannes Köster
'''

import re, sys, os, traceback, logging, glob
from multiprocessing import Pool, Event
from collections import defaultdict

from snakemake.rules import Rule
from snakemake.exceptions import MissingOutputException, MissingInputException, AmbiguousRuleException, CyclicGraphException, MissingRuleException, RuleException, CreateRuleException
from snakemake.utils import shell
from snakemake.jobs import Job, protected

class Workflow:
	def __init__(self):
		"""
		Create the controller.
		"""
		self.__rules = dict()
		self.__last = None
		self.__first = None
		self.__workdir_set = False
		self._jobs_finished = Event()
		self._virgin_globals = None
		self._runtimes = defaultdict(list)		
	
	def report_runtime(self, rule, runtime):
		self._runtimes[rule].append(runtime)
		
	def get_runtimes(self):
		for rule, runtimes in self._runtimes.items():
			s = sum(runtimes)
			yield rule, min(runtimes), max(runtimes), s, s / len(runtimes)
		
	def clear(self):
		self.__rules.clear()
		self.__last = None
		self.__first = None
		self.__workdir_set = False
		self._jobs_finished = Event()
		for k in list(globals().keys()):
			if k not in self._virgin_globals:
				del globals()[k]

	def setup_pool(self, jobs):
		self.__pool = Pool(processes=jobs)
		
	def set_jobs_finished(self, job = None):
		self._jobs_finished.set()
		
	def get_snakefile_globals(self):
		return self.__snakefile_globals
	
	def get_pool(self):
		"""
		Return the current thread pool.
		"""
		return self.__pool
	
	def add_rule(self, name):
		"""
		Add a rule.
		"""
		if self.is_rule(name):
			raise CreateRuleException("The name {} is already used by another rule".format(name))
		if "__" + name in globals():
			raise CreateRuleException("The name __{} is already used by a variable.".format(name))
		rule = Rule(name, self)
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

	def get_rule(self, name):
		"""
		Get rule by name.
		
		Arguments
		name -- the name of the rule
		"""
		return self.__rules[name]

	def last_rule(self):
		"""
		Return the last rule.
		"""
		return self.__last

	def run_first_rule(self, dryrun = False, forcethis = False, forceall = False):
		"""
		Apply the rule defined first.
		"""
		self.run_rule(self.__first, dryrun = dryrun, forcethis = forcethis, forceall = forceall)
		
	def run_rule(self, name, dryrun = False, forcethis = False, forceall = False):
		"""
		Apply a rule.
		
		Arguments
		name -- the name of the rule to apply
		"""
		rule = self.__rules[name]
		self._run(rule, forcethis = forcethis, forceall = forceall, dryrun = dryrun)
			
	def produce_file(self, file, dryrun = False, forcethis = False, forceall = False):
		"""
		Apply a rule such that the requested file is produced.
		
		Arguments
		file -- the path of the file to produce
		"""
		producer = None
		missing_input_ex = []
		for rule in self.__rules.values():
			if rule.is_producer(file):
				try:
					rule.run(file, jobs=dict(), forceall = forceall, dryrun = True, quiet = True, visited = set())
					if producer:
						raise AmbiguousRuleException("Ambiguous rules: {} and {}".format(producer, rule))
					producer = rule
				except MissingInputException as ex:
					missing_input_ex.append(ex)
		
		if not producer:
			if missing_input_ex:
				raise MissingInputException(include = missing_input_ex)
			raise MissingRuleException(file)
		self._run(producer, file, forcethis = forcethis, forceall = forceall, dryrun = dryrun)
	
	def _run(self, rule, requested_output = None, dryrun = False, forcethis = False, forceall = False):
		job = rule.run(requested_output, jobs=dict(), forcethis = forcethis, forceall = forceall, dryrun = dryrun, visited = set())
		job.run(callback = self.set_jobs_finished)
		self._jobs_finished.wait()

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

	def execdsl(self, compiled_dsl_code, rowmap):
		"""
		Execute a piece of compiled snakemake DSL.
		"""
		self.rowmap = rowmap
		exec(compiled_dsl_code, globals())

	def set_workdir(self, workdir):
		if not self.__workdir_set:
			if not os.path.exists(workdir):
				os.makedirs(workdir)
			os.chdir(workdir)
			self.__workdir_set = True

workflow = Workflow()

def _set_workdir(path):
	workflow.set_workdir(path)

def _add_rule(name):
	workflow.add_rule(name)

def _set_input(*paths, **kwpaths):
	workflow.last_rule().set_input(*paths, **kwpaths)

def _set_output(*paths, **kwpaths):
	workflow.last_rule().set_output(*paths, **kwpaths)

def _set_message(message):
	workflow.last_rule().set_message(message)

workflow._virgin_globals = dict(globals())
