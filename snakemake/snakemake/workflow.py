# -*- coding: utf-8 -*-

import re, os
from collections import OrderedDict
from itertools import chain

from snakemake.logging import logger
from snakemake.rules import Rule, Ruleorder
from snakemake.exceptions import MissingOutputException, MissingInputException, \
	AmbiguousRuleException, CyclicGraphException, MissingRuleException, \
	RuleException, CreateRuleException, ProtectedOutputException, \
	UnknownRuleException, NoRulesException, print_exception
from snakemake.shell import shell
from snakemake.dag import DAG
from snakemake.scheduler import JobScheduler
from snakemake.parser import compile_to_python
from snakemake.io import protected, temp, temporary, expand, dynamic, IOFile


__author__ = "Johannes Köster"

class Workflow:
	def __init__(self, snakemakepath = None):
		"""
		Create the controller.
		"""
		self._rules = OrderedDict()
		self.first_rule = None
		self._workdir = None
		self._ruleorder = Ruleorder()
		self.linemaps = dict()
		self.rule_count = 0
		self.snakemakepath = snakemakepath
		self.globals = globals()
	
	@property
	def rules(self):
		return self._rules.values()
	
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
		if not self.first_rule:
			self.first_rule = rule.name
		return name
			
	def is_rule(self, name):
		"""
		Return True if name is the name of a rule.
		
		Arguments
		name -- a name
		"""
		return name in self._rules

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
	
	def list_rules(details = True, log = logger.info):
		log("Available rules:")
		for rule in workflow.rules: 
			log(rule.name)
			if details:
				if rule.docstring:
					for line in rule.docstring.split("\n"):
						log("\t" + line)

	def execute(self, targets = None, dryrun = False,  touch = False, cores = 1,
	              forcetargets = False, forceall = False, forcerules = None, prioritytargets = None,
	              quiet = False, keepgoing = False, printshellcmds = False, printreason = False, printdag = False,
	              cluster = None,  ignore_ambiguity = False, workdir = None, stats = None):
		if workdir is None:
			workdir = os.getcwd() if self._workdir is None else self._workdir
		os.chdir(workdir)

		if not targets:
			targets = [self.first_rule]
		if not prioritytargets:
			prioritytargets = list()
		
		targetrules, targetfiles, priorityrules, priorityfiles = set(), set(), set(), set()
		for target in targets:
			if self.is_rule(target):
				rule = self._rules[target]
				targetrules.add(rule)
			else:
				file = os.path.relpath(target)
				targetfiles.add(file)
				
		for target in prioritytargets:
			if self.is_rule(target):
				rule = self._rules[target]
				priorityrules.add(rule)
				if not rule.has_wildcards():
					targetrules.add(rule)
			else:
				file = os.path.relpath(target)
				priorityfiles.add(file)
				targetfiles.add(file)
		
		try:
			forcerules_ = list()
			if forcerules:
				for r in forcerules:
					forcerules_.append(self._rules[r])
		except KeyError as ex:
			logger.critical("Rule {} is not available.".format(r))
			self.list_rules()
			return False
		
		dag = DAG(self, targetfiles=targetfiles, targetrules=targetrules, forceall=forceall, forcetargets=forcetargets, 
		          forcerules=forcerules_, priorityfiles=priorityfiles, priorityrules=priorityrules, ignore_ambiguity=ignore_ambiguity)
		try:
			dag.init()
		except RuleException as ex:
			# TODO think about printing a DAG on request here
			print_exception(ex, self.linemaps)
			return False
		except (Exception, BaseException) as ex:
			print_exception(ex, self.linemaps)
			return False
		
		if printdag:
			print(dag)
			return True
		
		scheduler = JobScheduler(self, dag, cores, dryrun=dryrun, touch=touch, cluster=cluster,
		                         quiet=quiet, keepgoing=keepgoing,
		                         printreason=printreason, printshellcmds=printshellcmds)
		success = scheduler.schedule()
		
		if success:
			if not dryrun and stats:
				scheduler.stats.to_csv(stats)
		else:
			logger.critical("Exiting because a job execution failed. Look above for error message")
			return False
		return True

	def include(self, snakefile, workdir = None, overwrite_first_rule = False):
		"""
		Include a snakefile.
		"""
		global workflow
		workflow = self
		first_rule = self.first_rule
		if workdir:
			os.chdir(workdir)
		code, linemap, rule_count = compile_to_python(snakefile, rule_count = self.rule_count)
		self.rule_count += rule_count
		self.linemaps[snakefile] = linemap
		exec(compile(code, snakefile, "exec"), self.globals)
		if not overwrite_first_rule:
			self.first_rule = first_rule

	def workdir(self, workdir):
		if self._workdir is None:
			if not os.path.exists(workdir):
				os.makedirs(workdir)
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
				rule.threads = ruleinfo.threads
			if ruleinfo.priority:
				rule.priority = ruleinfo.priority
			if ruleinfo.log:
				rule.log = ruleinfo.log
			if ruleinfo.message:
				rule.message = ruleinfo.message
			rule.docstring = ruleinfo.docstring
			rule.run_func = ruleinfo.func
			rule.shellcmd = ruleinfo.shellcmd
			return ruleinfo.func
		return decorate

	def docstring(self, string):
		def decorate(ruleinfo):
			ruleinfo.docstring = string
			return ruleinfo
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

	def priority(self, priority):
		def decorate(ruleinfo):
			ruleinfo.priority = priority
			return ruleinfo
		return decorate

	def log(self, log):
		def decorate(ruleinfo):
			ruleinfo.log = log
			return ruleinfo
		return decorate

	def shellcmd(self, cmd):
		def decorate(ruleinfo):
			ruleinfo.shellcmd = cmd
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
		self.shellcmd = None
		self.input = None
		self.output = None
		self.message = None
		self.threads = None
		self.priority = None
		self.log = None
		self.docstring = None
