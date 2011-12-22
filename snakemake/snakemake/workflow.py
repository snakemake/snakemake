# -*- coding: utf-8 -*-

'''
Created on 13.11.2011

@author: Johannes Köster
'''

import re, os, logging, subprocess, glob, inspect
from multiprocessing import Pool
from collections import defaultdict



# Global functions
if "SHELL" in os.environ:
	def _shell(cmd):
		subprocess.check_call(cmd, shell=True, executable = os.environ["SHELL"])
else:
	def _shell(cmd):
		subprocess.check_call(cmd, shell=True)

def shell(cmd, *args, **kwargs):
	variables = dict(globals())
	# add local variables from calling rule/function
	variables.update(inspect.currentframe().f_back.f_locals)
	variables.update(kwargs)
	_shell(cmd.format(*args, **variables))

class RuleException(Exception):
	pass

def run_wrapper(run, rulename, ruledesc, input, output, wildcards):
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
		raise RuleException(": ".join(type(ex).__name__,str(ex)))
	for o in output:
		if not os.path.exists(o):
			raise RuleException("Output file {} not produced by rule {}.".format(o, rulename))

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
	
	def _expand_wildcards(self, requested_output):
		if requested_output:
			wildcards = self.get_wildcards(requested_output[0])
		else:
			return tuple(self.input), tuple(self.output), dict()

		try:
			input = tuple(i.format(**wildcards) for i in self.input)
			output = tuple(o.format(**wildcards) for o in self.output)
			return input, output, wildcards
		except KeyError:
			raise RuleException("Could not resolve wildcard in rule {}: {}".format(self.name, i))

	def _to_visit(self, input, forceall = False):
		if forceall:
			missing_input = input
		else:
			missing_input = tuple(i for i in input if not os.path.exists(i))
		rules = workflow.get_rules()

		producer = defaultdict(list)
		noproducer = set(missing_input)
		for rule in rules:
			if rule != self:
				for i in missing_input:
					if rule.is_producer(i):
						producer[rule].append(i)
						noproducer.remove(i)

		if noproducer:
			raise RuleException("Missing input files in rule {}:\n{}.".format(self.name, ", ".join(noproducer)))

		tovisit = dict()
		for rule, files in producer.items():
			for request_output in rule.partition_output(files):
				tovisit[rule] = request_output
		return tovisit
		
	
	def check_dag(self, requested_output = [], forceall = False, visited = set()):
		visited.add(self)
		nodes = 1

		input, output, _ = self._expand_wildcards(requested_output)

		tovisit = self._to_visit(input, forceall = forceall)
		
		input_provider = dict()
		for rule, files in tovisit.items():
			for i in files:
				if i in input_provider:
					raise RuleException("Ambiguous rules: {} and {}".format(rule.name, input_provider[i]))
				elif rule in visited:
					raise RuleException("Circular dependency between {} and {}".format(self.name, rule.name))
				else:
					input_provider[i] = rule

		for rule in set(input_provider.values()):
			nodes += rule.check_dag(tovisit[rule], visited = set(visited))
		return nodes

	def run(self, requested_output = [], jobs = dict(), forcethis = False, forceall = False):
		input, output, wildcards = self._expand_wildcards(requested_output)
		tovisit = self._to_visit(input, forceall = forceall)

		todo = []
		for rule, files in tovisit.items():
			todo.append(rule.run(files, jobs, forceall = forceall))
		for job in todo:
			job.get()

		if forcethis or forceall or self._need_run(input, output, jobs):
			job = workflow.get_pool().apply_async(
					run_wrapper, 
					[self._get_run(), self.name, self._get_message(input, output, wildcards), input, output, wildcards])
			jobs[output] = job
			return job

	def dryrun(self, requested_output = [], jobs = set(), forcethis = False, forceall = False):
		input, output, wildcards = self._expand_wildcards(requested_output)
		tovisit = self._to_visit(input, forceall = forceall)

		for rule, files in tovisit.items():
			rule.dryrun(files, jobs, forceall = forceall)

		if forcethis or forceall or self._need_run(input, output, jobs):
			print(self._get_message(input, output, wildcards))
			jobs.add(output)

	def check(self):
		if self.output and not self.has_run():
			raise RuleException("Rule {} defines output but does not have a \"run\" definition.".format(self.name))

	def _need_run(self, input, output, jobs):
		if not self.has_run():
			return False
		if output in jobs:
			return False
		for o in output:
			if not os.path.exists(o): return True
		mintime = min(map(lambda f: os.stat(f).st_mtime, output))
		for i in input:
			if os.stat(i).st_mtime >= mintime: return True
		return False

	def _get_run(self):
		return globals()[self.name]

	def has_run(self):
		return self.name in globals()

	def _get_message(self, input, output, wildcards):
		return "rule {name}:\n\tinput: {input}\n\toutput: {output}\n".format(
			name=self.name, input=", ".join(input), output=", ".join(output))

	def is_parent(self, rule):
		return self in rule.parents.values()

	def is_producer(self, requested_output):
		"""
		Returns True if this rule is a producer of the requested output.
		"""
		for o in self.regex_output:
			match = re.match(o, requested_output)
			if match and len(match.group()) == len(requested_output):
				return True
		return False
			
	def _wildcards_to_str(self, wildcards):
		if wildcards:
			return "Wildcards:\n" + "\n".join(": ".join(i) for i in wildcards.items())
		return ""

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

class Workflow:

	def __init__(self):
		"""
		Create the controller.
		"""
		self.__rules = dict()
		self.__last = None
		self.__first = None
		self.__workdir_set = False

	def setup_pool(self, jobs):
		self.__pool = Pool(processes=jobs)
	
	def get_pool(self):
		"""
		Return the current thread pool.
		"""
		return self.__pool
	
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

	def run_first_rule(self, dryrun = False, forcethis = False, forceall = False):
		"""
		Apply the rule defined first.
		"""
		self.__first.check_dag()
		if dryrun:
			self.__first.dryrun(forcethis = forcethis, forceall = forceall)
		else:
			job = self.__first.run(forcethis = forcethis, forceall = forceall)
			if job: job.get()		
		
	def run_rule(self, name, dryrun = False, forcethis = False, forceall = False):
		"""
		Apply a rule.
		
		Arguments
		name -- the name of the rule to apply
		"""
		rule = self.__rules[name]
		if dryrun:
			rule.dryrun(forcethis = forcethis, forceall = forceall)
		else:
			job = rule.run(forcethis = forcethis, forceall = forceall)
			if job: job.get()

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

workflow = Workflow()

def _set_workdir(path):
	workflow.set_workdir(path)

def _add_rule(name):
	workflow.add_rule(Rule(name))

def _set_input(paths):
	workflow.last_rule().add_input(paths)

def _set_output(paths):
	workflow.last_rule().add_output(paths)
