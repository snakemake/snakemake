# -*- coding: utf-8 -*-

'''
Created on 13.11.2011

@author: Johannes Köster
'''

import re, os, logging, subprocess, glob
from multiprocessing import Pool
from collections import defaultdict


# Global functions
def shell(cmd, *args, **kwargs):
	cmd = cmd.format(*args, **kwargs)
	if "SHELL" in os.environ:
		subprocess.check_call(cmd, shell=True, executable = os.environ["SHELL"])
	else:
		subprocess.check_call(cmd, shell=True)

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
		raise RuleException(str(ex))

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
		self.wildcard_names = set()
		self.jobs = dict()

	def __to_regex(self, output):
		"""
		Convert a filepath containing wildcards to a regular expression.
		
		Arguments
		output -- the filepath
		"""
		output = re.sub("\.", "\.", output)
		return re.sub('\{(?P<name>\w+?)\}', lambda match: '(?P<{}>.+)'.format(match.group('name')), output)

	def _get_wildcard_names(self, output):
		return set(match.group('name') for match in re.finditer("\{(?P<name>\w+?)\}", output))

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
				wildcards = self._get_wildcard_names(item)
				if self.output:
					if self.wildcard_names != wildcards:
						raise RuleException("Not all output files of rule {} contain the same wildcards. ".format(self.name))
				else:
					self.wildcard_names = wildcards
				self.output.append(item)
				self.regex_output.append(self.__to_regex(item))

	def is_parent(self, rule):
		return self in rule.parents.values()

	def setup_parents(self, wildcards = {}, requested_output = []):
		"""
		Setup the DAG by finding parent rules that create files needed as input for this rule
		"""
		if not wildcards and requested_output:
			wildcards = self.get_wildcards(requested_output[0])

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
		
		for rule, files in products.items():
			for partition in rule.partition_output(files):
				rule.setup_parents(dict(), partition)

	def is_producer(self, requested_output):
		"""
		Returns True if this rule is a producer of the requested output.
		"""
		for o in self.regex_output:
			match = re.match(o, requested_output)
			if match and len(match.group()) == len(requested_output):
				return True
		return False

	def check_input(self, files, wildcards):
		"""
		Check if all input files present.
		"""
		notpresent = [f for f in files if not os.path.exists(f)]
		if notpresent:
			raise RuleException("Missing input files for rule {}: {}\n{}".format(
				self.name, ", ".join(notpresent), self._wildcards_to_str(wildcards)))
			
	def _wildcards_to_str(self, wildcards):
		if wildcards:
			return "Wildcards:\n" + "\n".join(": ".join(i) for i in wildcards.items())
		return ""

	def expand_input(self, input, flat = True):
		"""
		Expand unix wildcards in input files.
		"""
		expand = lambda input, expanded: input.append(expanded)
		if flat:
			expand = lambda input, expanded: input.extend(expanded)

		input = list(input)
		for i in input:
			expanded = glob.glob(i)
			if len(expanded) == 0:
				expanded = i
			elif len(expanded) == 1:
				expanded = expanded[0]
			expand(input, expanded)
		return input

	def get_wildcards(self, requested_output):
		"""
		Update the given wildcard dictionary by matching regular expression output files to the requested concrete ones.
		
		Arguments
		wildcards -- a dictionary of wildcards
		requested_output -- a concrete filepath
		"""
		bestmatchlen = 0
		bestmatch = None
		for o in self.regex_output:
			match = re.match(o, requested_output)
			if match and len(match.group()) == len(requested_output):
				l = self.get_wildcard_len(match.groupdict())
				if not bestmatch or bestmatchlen > l:
					bestmatch = match.groupdict()
					bestmatchlen = l
		return bestmatch
		
	
	def get_wildcard_len(self, wildcards):
		return sum(map(len, wildcards.values()))

	def partition_output(self, requested_outputs):
		partition = defaultdict(list)
		for r in requested_outputs:
			wc = frozenset(self.get_wildcards(r).items())
			partition[wc].append(r)
		return partition.values()
			
	def apply_rule(self, wildcards = {}, requested_output = [], dryrun = False, force =False):
		"""
		Apply the rule
		
		Arguments
		wildcards -- a dictionary of wildcards
		requested_output -- the requested concrete output file 
		"""
		if not wildcards and requested_output:
			wildcards = self.get_wildcards(requested_output[0]) # wildcards can be determined with only one output since each has to use the same		

		output = [o.format(**wildcards) for o in self.output]
		input = [i.format(**wildcards) for i in self.input]

		products = defaultdict(list)
		notproduced = []
		for i in input:
			if i in self.parents:
				products[self.parents[i]].append(i)
			else:
				notproduced.append(i)

		jobs = []
		for rule, files in products.items():
			for partition in rule.partition_output(files):
				jobs.append((rule, rule.apply_rule(requested_output = partition, dryrun = dryrun, force = force)))
		Controller.get_instance().join_pool(jobs = jobs)

		# all inputs have to be present after finishing parent jobs
		if dryrun:
			self.check_input(notproduced, wildcards)
		else:
			self.check_input(input, wildcards)
			
		if len(output) > 0 and Controller.get_instance().is_produced(output):
			# if output is already produced, only recalculate if input is newer.
			time = min(map(lambda f: os.stat(f).st_mtime, output))
			if not force and Controller.get_instance().is_produced(input) and not Controller.get_instance().is_newer(input, time):
				return

		if self.name in globals():
			_output = tuple(output)
			if dryrun:
				if not _output in self.jobs:
					# print job if not yet printed
					self.print_job(input, output)
					self.jobs[_output] = None
			else:
				if _output in self.jobs:
					# job already started
					return self.jobs[_output]
					
				self.print_job(input, output)
				# if there is a run body, run it asyncronously
				job = Controller.get_instance().get_pool().apply_async(run_wrapper, [globals()[self.name], input, output, wildcards])
				self.jobs[_output] = job
				return job
	
	def print_job(self, input, output):
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
