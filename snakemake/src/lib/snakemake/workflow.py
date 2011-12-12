# -*- coding: utf-8 -*-

'''
Created on 13.11.2011

@author: Johannes Köster
'''

import re, os, logging
from multiprocessing import Pool
from collections import defaultdict

class RuleException(Exception):
	pass

def run_wrapper(run, input, output, wildcards):
	"""
	Wrapper around the run method that handles directory creation and output file deletion on error.
	
	Arguments
	run -- the run method
	input -- list of input files
	output -- list of output files
	wildcards -- so far processed wildcards
	"""
	for o in output:
		dir = os.path.dirname(o)
		if len(dir) > 0 and not os.path.exists(dir):
			os.makedirs(dir)
	try:
		run(input, output, wildcards)
	except (Exception, BaseException) as ex:
		# Remove produced output on exception
		for o in output:
			if os.path.isdir(o): os.rmdir(o)
			elif os.path.exists(o): os.remove(o)
		raise ex

class Rule:
	def __init__(self, name):
		"""
		Create a rule
		
		Arguments
		name -- the name of the rule
		"""
		self.name = name
		self.input = []
		self.output = []
		self.regex_output = []
		self.parents = dict()

	def __to_regex(self, output):
		"""
		Convert a filepath containing wildcards to a regular expression.
		
		Arguments
		output -- the filepath
		"""
		return re.sub('\{(?P<name>.+?)\}', lambda match: '(?P<{}>.+?)'.format(match.group('name')), output)

	def add_input(self, input):
		"""
		Add a list of input files. Recursive lists are flattened.
		
		Arguments
		input -- the list of input files
		"""
		for item in input:
			if isinstance(item, list): self.add_input(item)
			else: self.input.append(item)

	def add_output(self, output):
		"""
		Add a list of output files. Recursive lists are flattened.
		
		Arguments
		output -- the list of output files
		"""
		for item in output:
			if isinstance(item, list): self.add_output(item)
			else:
				self.output.append(item)
				self.regex_output.append(self.__to_regex(item))

	def is_parent(self, rule):
		return self in rule.parents.values()

	def setup_parents(self, wildcards = {}, requested_output = []):
		"""
		Setup the DAG by finding parent rules that create files needed as input for this rule
		"""
		for o in requested_output:
			self.update_wildcards(wildcards, o)

		products = defaultdict(list)
		for i in self.input:
			try:
				i = i.format(**wildcards)
			except KeyError:
				raise RuleException("Could not resolve wildcard in rule {}: {}".format(self.name, i))
			found = None
			for rule in Controller.get_instance().get_rules():
				if rule != self and rule.is_producer(i):
					if self.is_parent(rule):
						raise IOError("Circular dependency between rules: {} and {}".format(rule.name, self.name))
					if found:
						raise IOError("Ambiguous rules: {} and {}".format(rule.name, found))
					self.parents[i] = rule
					products[rule].append(i)
					found = rule.name

		for rule in products:
			rule.setup_parents(dict(wildcards), products[rule])

	def is_producer(self, requested_output):
		"""
		Returns True if this rule is a producer of the requested output.
		"""
		for o in self.regex_output:
			match = re.match(o, requested_output)
			if match: return True
		return False

	def update_wildcards(self, wildcards, requested_output):
		"""
		Update the given wildcard dictionary by matching regular expression output files to the requested concrete ones.
		
		Arguments
		wildcards -- a dictionary of wildcards
		requested_output -- a concrete filepath
		"""
		wildcards = dict(wildcards)
		for o in self.regex_output:
			match = re.match(o, requested_output)
			if match:
				wildcards.update(match.groupdict())
				return wildcards
		return wildcards

	def apply_rule(self, wildcards = {}, requested_output = [], dryrun = False, force =False):
		"""
		Apply the rule
		
		Arguments
		wildcards -- a dictionary of wildcards
		requested_output -- the requested concrete output file 
		"""
		for o in requested_output:
			wildcards = self.update_wildcards(wildcards, o)

		output = [o.format(**wildcards) for o in self.output]
		input = [i.format(**wildcards) for i in self.input]

		products = defaultdict(list)
		for i in input:
			if i in self.parents:
				products[self.parents[i]].append(i)
		jobs = []
		for rule, files in products.items():
			jobs.append((rule, rule.apply_rule(wildcards, files, dryrun = dryrun, force = force)))
		Controller.get_instance().join_pool(jobs = jobs)

		# all inputs have to be present after finishing parent jobs
		if not dryrun and not Controller.get_instance().is_produced(input):
			raise RuleException("Could not execute rule {}: not all input files present.".format(self.name))
			
		if len(output) > 0 and Controller.get_instance().is_produced(output):
			# if output is already produced, only recalculate if input is newer.
			time = min(map(lambda f: os.stat(f).st_mtime, output))
			if not force and Controller.get_instance().is_produced(input) and not Controller.get_instance().is_newer(input, time):
				return

		self.print_rule(input, output)
		if not dryrun and self.name in globals():
			# if there is a run body, run it asyncronously
			return Controller.get_instance().get_pool().apply_async(run_wrapper, [globals()[self.name], input, output, wildcards])
	
	def print_rule(self, input, output):
		 logging.info("rule {name}:\n\tinput: {input}\n\toutput: {output}\n".format(name=self.name, input=", ".join(input), output=", ".join(output)))

class Controller:
	instance = None
	jobs = 1

	@classmethod
	def get_instance(cls):
		"""
		Return the singleton instance of the controller.
		"""
		if cls.instance:
			return cls.instance
		cls.instance = cls()
		return cls.instance

	def __init__(self):
		"""
		Create the controller.
		"""
		self.__rules = dict()
		self.__last = None
		self.__first = None
		self.__workdir_set = False

	def setup_pool(self):
		self.__pool = Pool(processes=Controller.jobs)
	
	def get_pool(self):
		"""
		Return the current thread pool.
		"""
		return self.__pool


	def join_pool(self, jobs = None, rule = None, result = None):
		"""
		Join all threads in pool together.
		"""
		if jobs:
			for rule, result in jobs:
				if result:
					try: result.get() # reraise eventual exceptions
					except (Exception, BaseException) as ex:
						raise RuleException("Error: Could not execute rule {}: {}".format(rule.name, str(ex)))
		elif rule and result:
			try: result.get() # reraise eventual exceptions
			except (Exception, BaseException) as ex:
				raise RuleException("Error: Could not execute rule {}: {}".format(rule.name, str(ex)))
	
	def add_rule(self, rule):
		"""
		Add a rule.
		"""
		self.__rules[rule.name] = rule
		self.__last = rule
		if not self.__first:
			self.__first = rule
			
	def is_rule(self, name):
		"""
		Return True if name is the name of a rule.
		
		Arguments
		name -- a name
		"""
		return name in self.__rules

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

	def apply_first_rule(self, dryrun = False, force = False):
		"""
		Apply the rule defined first.
		"""
		self.__first.setup_parents()
		self.setup_pool()
		self.join_pool(rule = self.__first, result = self.__first.apply_rule(dryrun = dryrun, force = force))
		
		
	def apply_rule(self, name, dryrun = False, force = False):
		"""
		Apply a rule.
		
		Arguments
		name -- the name of the rule to apply
		"""
		self.__rules[name].setup_parents()
		self.setup_pool()
		self.join_pool(rule = self.__rules[name], result = self.__rules[name].apply_rule(dryrun = dryrun, force = force))

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

	def execdsl(self, compiled_dsl_code):
		"""
		Execute a piece of compiled snakemake DSL.
		"""
		exec(compiled_dsl_code, globals())

	def set_workdir(self, workdir):
		if not self.__workdir_set:
			if not os.path.exists(workdir):
				os.makedirs(workdir)
			os.chdir(workdir)
			self.__workdir_set = True

def _set_workdir(path):
	Controller.get_instance().set_workdir(path)

def _add_rule(name):
	Controller.get_instance().add_rule(Rule(name))

def _set_input(paths):
	Controller.get_instance().last_rule().add_input(paths)

def _set_output(paths):
	Controller.get_instance().last_rule().add_output(paths)
