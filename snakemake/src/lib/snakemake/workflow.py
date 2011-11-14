# -*- coding: utf-8 -*-

import re
import os
from multiprocessing import Pool

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
	except:
		for o in output:
			if os.path.isdir(o): os.rmdir(o)
			else: os.remove(o)

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

	def setup_parents(self):
		"""
		Setup the DAG by finding parent rules that create files needed as input for this rule
		"""
		for i in self.input:
			found = None
			for rule in Controller.get_instance().get_rules():
				if rule.is_producer(i):
					if found:
						raise IOError("Ambiguous rules: {} and {}".format(rule.name, found))
					self.parents[i] = rule
					found = rule.name

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

	def apply_rule(self, wildcards = {}, requested_output = None):
		"""
		Apply the rule
		
		Arguments
		wildcards -- a dictionary of wildcards
		requested_output -- the requested concrete output file 
		"""
		if requested_output:
			wildcards = self.update_wildcards(wildcards, requested_output)

		output = [o.format(**wildcards) for o in self.output]
		
		if Controller.get_instance().is_produced(output):
			return

		input = [i.format(**wildcards) for i in self.input]

		for i in range(len(self.input)):
			if self.input[i] in self.parents:
				self.parents[self.input[i]].apply_rule(wildcards, input[i])
		Controller.get_instance().join_pool()
		
		# all inputs have to be present after finishing parent jobs
		if not Controller.get_instance().is_produced(input):
			exit(1)

		Controller.get_instance().get_pool().apply(run_wrapper, [globals()[self.name], input, output, wildcards]) 

class Controller:
	instance = None
	processes = 1

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
		self.__pool = Pool(processes=Controller.processes)
	
	def get_pool(self):
		"""
		Return the current thread pool.
		"""
		return self.__pool
	
	def join_pool(self):
		"""
		Join all threads in pool together.
		"""
		self.__pool.close()
		self.__pool.join()
		self.__pool = Pool(processes=Controller.processes)
	
	def add_rule(self, rule):
		"""
		Add a rule.
		"""
		self.__rules[rule.name] = rule
		self.__last = rule
		if not self.__first:
			self.__first = rule

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

	def apply_first_rule(self):
		"""
		Apply the rule defined first.
		"""
		self.__first.apply_rule()
		
	def apply_rule(self, name):
		"""
		Apply a rule.
		
		Arguments
		name -- the name of the rule to apply
		"""
		self.__rules[name].apply_rule()

	def get_rules(self):
		"""
		Get the list of rules.
		"""
		return self.__rules.values()

	def setup_dag(self):
		"""
		Setup the DAG.
		"""
		for rule in self.get_rules():
			rule.setup_parents()

	def is_produced(self, output):
		"""
		Return True if file is already produced.
		"""
		for o in output:
			if not os.path.exists(o): return False
		return True

	def execdsl(self, compiled_dsl_code):
		"""
		Execute a piece of compiled snakemake DSL.
		"""
		exec(compiled_dsl_code, globals())