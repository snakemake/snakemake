# -*- coding: utf-8 -*-

import re, sys, os, traceback, glob, signal
from multiprocessing import Event
from collections import defaultdict, OrderedDict
from itertools import chain
from functools import lru_cache
from tempfile import TemporaryFile

from snakemake.logging import logger
from snakemake.rules import Rule, Ruleorder
from snakemake.exceptions import MissingOutputException, MissingInputException, \
	AmbiguousRuleException, CyclicGraphException, MissingRuleException, \
	RuleException, CreateRuleException, ProtectedOutputException, \
	UnknownRuleException, NoRulesException
from snakemake.shell import shell, format
from snakemake.jobs import Job, KnapsackJobScheduler, ClusterJobScheduler, print_job_dag
from snakemake.parser import compile_to_python
from snakemake.io import protected, temp, temporary, expand, dynamic, IOFile


__author__ = "Johannes Köster"

class Jobcounter:
	def __init__(self, count):
		self.count = count
		self._done = 0
	
	def done(self):
		self._done += 1
	
	def __str__(self):
		return "{} of {} steps ({}%) done".format(self._done, self.count, int(self._done / self.count * 100))

class Workflow:
	def __init__(self):
		"""
		Create the controller.
		"""
		self._rules = OrderedDict()
		self._first = None
		self._workdir = None
		self._runtimes = defaultdict(list)
		self._ruleorder = Ruleorder()
		self.cores = 1
		self.rowmaps = dict()
		self.jobcounter = None
		self.rule_count = 0
		self.errors = False
	
	def report_runtime(self, rule, runtime):
		self._runtimes[rule].append(runtime)

	def get_ruleorder(self):
		return self._ruleorder
		
	def get_runtimes(self):
		for rule, runtimes in self._runtimes.items():
			s = sum(runtimes)
			yield rule, min(runtimes), max(runtimes), s, s / len(runtimes)

	def get_rule_count(self):
		return len(self._rules)
	
	def add_rule(self, name = None, lineno = None, snakefile = None):
		"""
		Add a rule.
		"""
		if name == None:
			name = str(len(self._rules))
		if self.is_rule(name):
			raise CreateRuleException(
				"The name {} is already used by another rule".format(name))
		rule = Rule(name, self, lineno = lineno, snakefile = snakefile)
		self._rules[rule.name] = rule
		if not self._first:
			self._first = rule.name
		return name
			
	def is_rule(self, name):
		"""
		Return True if name is the name of a rule.
		
		Arguments
		name -- a name
		"""
		return name in self._rules
	
	def get_producers(self, files, exclude = None):
		checked = set()
		for f in files:
			if not f in checked:
				checked.add(f)
				for item in self._get_producers(f, exclude=exclude):
					yield item

	@lru_cache()
	def _get_producers(self, file, exclude = None):
		producers = []
		for rule in self.get_rules():
			if rule != exclude:
				if rule.is_producer(file):
					producers.append((rule, file))
		return producers

	def get_rule(self, name):
		"""
		Get rule by name.
		
		Arguments
		name -- the name of the rule
		"""
		if not self._rules:
			raise NoRulesException()
		if not name in self._rules:
			raise UnknownRuleException(name)
		return self._rules[name]

	def run_first_rule(self, dryrun = False, touch = False, 
		forcethis = False, forceall = False, give_reason = False, 
		cluster = None, dag = False, ignore_ambiguity = False):
		"""
		Apply the rule defined first.
		"""
		first = self._first
		if not first:
			for key, value in self._rules.items():
				first = key
				break
		return self._run([(self.get_rule(first), None)], 
			dryrun = dryrun, touch = touch, forcethis = forcethis, 
			forceall = forceall, give_reason = give_reason, 
			cluster = cluster, dag = dag, ignore_ambiguity = ignore_ambiguity)
			
	def get_file_producers(self, files, dryrun = False, 
		forcethis = False, forceall = False, ignore_ambiguity = False):
		"""
		Return a dict of rules with requested files such that the requested files are produced.
		
		Arguments
		files -- the paths of the files to produce
		"""
		producers = dict()
		missing_input_ex = defaultdict(list)
		for rule, file in self.get_producers(files):
			try:
				job = rule.run(file, jobs=dict(), forceall = forceall, 
					dryrun = True, visited = set())
				if file in producers:
					if producers[file].rule > rule:
						continue
					elif producers[file].rule < rule:
						pass
					else:
						if ignore_ambiguity:
							logger.warning("Rules {rule1} and {} are ambigous for file {}, using {rule1}.".format(rule, file, rule1=produced[file].rule))
							continue
						raise AmbiguousRuleException(file, producers[file], job,
						                             lineno = rule.lineno,
						                             snakefile = rule.snakefile)
				producers[file] = job
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

		# clear eventually created IOFiles. especially needed to empty the count for needed files
		IOFile.clear()

		return [(job.rule, file) for file, job in producers.items()]

	def run_rules(self, targets, dryrun = False, touch = False, 
		forcethis = False, forceall = False, give_reason = False, 
		cluster = None, dag = False, ignore_ambiguity = False):
		ruletargets, filetargets = [], []
		for target in targets:
			if workflow.is_rule(target):
				ruletargets.append(target)
			else:
				filetargets.append(os.path.relpath(target))
		try:
			torun = self.get_file_producers(filetargets, forcethis = forcethis, 
				forceall = forceall, dryrun = dryrun, ignore_ambiguity = ignore_ambiguity) + \
				[(self.get_rule(name), None) for name in ruletargets]
		except AmbiguousRuleException as ex:
			if not dag:
				raise ex
			print_job_dag(chain(ex.job1.all_jobs(), ex.job2.all_jobs()))
			return
				
		return self._run(torun, dryrun = dryrun, touch = touch, 
			forcethis = forcethis, forceall = forceall, 
			give_reason = give_reason, cluster = cluster, dag = dag, 
			ignore_ambiguity = ignore_ambiguity)
	
	def _run(self, torun, dryrun = False, touch = False, forcethis = False, 
		forceall = False, give_reason = False, cluster = None, dag = False,
		ignore_ambiguity = False):
		jobs = dict()
		Job.count = 0
		
		root_jobs = set()
		try:
			for rule, requested_output in torun:
				root_jobs.add(rule.run(requested_output, jobs=jobs, forcethis = forcethis, 
					forceall = forceall, dryrun = dryrun, give_reason = give_reason, 
					touch = touch, visited = set(), 
					ignore_ambiguity = ignore_ambiguity))
		

			# collect all jobs
			all_jobs = set()
			for job in root_jobs:
				all_jobs.update(job.all_jobs())

			self.jobcounter = Jobcounter(sum(1 for job in all_jobs if not job.pseudo))
		except AmbiguousRuleException as ex:
			if not dag:
				raise ex
			else:
				all_jobs = chain(ex.job1.all_jobs(), ex.job2.all_jobs())
		
		if dag:
			print_job_dag(all_jobs)
			return
		if cluster:
			scheduler = ClusterJobScheduler(all_jobs, self, submitcmd = cluster)
		else:
			scheduler = KnapsackJobScheduler(all_jobs, self)
		success = scheduler.schedule()

		if not success:
			Job.cleanup_unfinished(all_jobs)
			logger.critical(
				"Exiting because a job execution failed. Look above for error message")
			return False
		return True

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
		return self._rules.values()

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
		global workflow
		workflow = self
		first_rule = self._first
		code, rowmap, rule_count = compile_to_python(snakefile, rule_count = self.rule_count)
		self.rule_count += rule_count
		self.rowmaps[snakefile] = rowmap
		exec(compile(code, snakefile, "exec"), globals())
		if not overwrite_first_rule:
			self._first = first_rule

	def workdir(self, workdir):
		if not self._workdir:
			if not os.path.exists(workdir):
				os.makedirs(workdir)
			os.chdir(workdir)
			self._workdir = workdir

	def ruleorder(self, *rulenames):
		self._ruleorder.add(*rulenames)

	def rule(self, name = None, lineno = None, snakefile = None):
		name = self.add_rule(name, lineno, snakefile)
		rule = self.get_rule(name)
		def decorate(ruleinfo):
			if ruleinfo.input:
				rule.set_input(*ruleinfo.input[0], **ruleinfo.input[1])
			if ruleinfo.output:
				rule.set_output(*ruleinfo.output[0], **ruleinfo.output[1])
			if ruleinfo.threads:
				rule.set_threads(ruleinfo.threads)
			if ruleinfo.message:
				rule.set_message(ruleinfo.message)
			rule.run_func = ruleinfo.func
			return ruleinfo.func
		return decorate


	def input(self, *paths, **kwpaths):
		def decorate(ruleinfo):
			ruleinfo.input = (paths, kwpaths)
			return ruleinfo
		return decorate

	def output(self, *paths, **kwpaths):
		def decorate(ruleinfo):
			ruleinfo.output = (paths, kwpaths)
			return ruleinfo
		return decorate

	def message(self, message):
		def decorate(ruleinfo):
			ruleinfo.message = message
			return ruleinfo
		return decorate

	def threads(self, threads):
		def decorate(ruleinfo):
			ruleinfo.threads = threads
			return ruleinfo
		return decorate

	def run(self, func):
		return RuleInfo(func)

	@staticmethod
	def _empty_decorator(f):
		return f


class RuleInfo:
	def __init__(self, func):
		self.func = func
		self.input = None
		self.output = None
		self.message = None
		self.threads = None
