import os, re
from operator import itemgetter
from collections import defaultdict
from snakemake.jobs import Job
from snakemake.exceptions import MissingInputException, AmbiguousRuleException, CyclicGraphException, RuleException

class Namedlist(list):
	def __init__(self, toclone = None, fromdict = None):
		super(Namedlist, self).__init__()
		self._names = dict()
		if toclone:
			self.extend(toclone)
			if isinstance(toclone, Namedlist):
				self.take_names(toclone.get_names())
		if fromdict:
			for key, item in fromdict.items():
				self.append(item)
				self.add_name(key)

	def add_name(self, name):
		self.set_name(name, len(self) - 1)
	
	def set_name(self, name, index):
		self._names[name] = index
		setattr(self, name, self[index])
			
	def get_names(self):
		for name, index in self._names.items():
			yield name, index
	
	def take_names(self, names):
		for name, index in names:
			self.set_name(name, index)
	
	def __hash__(self):
		return hash(tuple(self))

class Nameddict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

class Rule:
	def __init__(self, name, workflow):
		"""
		Create a rule
		
		Arguments
		name -- the name of the rule
		"""
		self.name = name
		self.message = None
		self.input = Namedlist()
		self.output = Namedlist()
		self.regex_output = []
		self.wildcard_names = set()
		self.workflow = workflow

	def _to_regex(self, output):
		"""
		Convert a filepath containing wildcards to a regular expression.
		
		Arguments
		output -- the filepath
		"""
		output = re.sub("\.", "\.", output)
		return re.sub('\{(?P<name>\w+?)\}', lambda match: '(?P<{}>.+)'.format(match.group('name')), output)

	def _get_wildcard_names(self, output):
		return set(match.group('name') for match in re.finditer("\{(?P<name>\w+?)\}", output))

	def has_wildcards(self):
		return bool(self.wildcard_names)

	def set_input(self, *input, **kwinput):
		"""
		Add a list of input files. Recursive lists are flattened.
		
		Arguments
		input -- the list of input files
		"""
		for item in input:
			self._set_inoutput_item(item, self.input)
		for name, item in kwinput.items():
			self._set_inoutput_item(item, self.input, name = name)

	def set_output(self, *output, **kwoutput):
		"""
		Add a list of output files. Recursive lists are flattened.
		
		Arguments
		output -- the list of output files
		"""
		for item in output:
			self._set_inoutput_item(item, self.output)
		for name, item in kwoutput.items():
			self._set_inoutput_item(item, self.output, name = name)
		
		for item in self.output:
			wildcards = self._get_wildcard_names(item)
			if self.wildcard_names:
				if self.wildcard_names != wildcards:
					raise SyntaxError("Not all output files of rule {} contain the same wildcards. ".format(self.name))
			else:
				self.wildcard_names = wildcards
			self.regex_output.append(self._to_regex(item))
	
	def _set_inoutput_item(self, item, inoutput, name=None):
		if isinstance(item, str):
			inoutput.append(item)
			if name:
				inoutput.add_name(name)
		else:
			for i in item:
				self._set_inoutput_item(i)

	def set_message(self, message):
		"""
		Set the message that is displayed when rule is executed.
		"""
		self.message = message
	
	def _expand_wildcards(self, requested_output):
		""" Expand wildcards depending on the requested output. """
		if requested_output:
			wildcards = self.get_wildcards(requested_output)
			missing_wildcards = set(wildcards.keys()) - self.wildcard_names 
		elif self.has_wildcards():				
			missing_wildcards = self.wildcard_names
		else:
			return Namedlist(self.input), Namedlist(self.output), dict()
		
		if missing_wildcards:
			raise RuleException("Could not resolve wildcards in rule {}:\n{}".format(self.name, "\n".join(self.wildcard_names)))

		try:
			input = Namedlist(i.format(**wildcards) for i in self.input)
			output = Namedlist(o.format(**wildcards) for o in self.output)
			input.take_names(self.input.get_names())
			output.take_names(self.output.get_names())
			return input, output, wildcards
		except KeyError as ex:
			# this can only happen if an input file contains an unresolved wildcard.
			raise SyntaxError("Wildcards in input file of rule {} do not appear in output files:\n{}".format(rule, str(ex)))
			

	def _get_missing_files(self, files):
		""" Return the tuple of files that are missing form the given ones. """
		return tuple(f for f in files if not os.path.exists(f))
	
	def _has_missing_files(self, files):
		""" Return True if any of the given files does not exist. """
		for f in files:
			if not os.path.exists(f):
				return True
		return False
	
	def _to_visit(self, input):
		for rule in self.workflow.get_rules():
			if rule != self:
				for i in input:
					if rule.is_producer(i):
						yield rule, i
	
	def run(self, requested_output = None, forceall = False, forcethis = False, jobs = dict(), dryrun = False, quiet = False, visited = set()):
		if (self, requested_output) in visited:
			raise CyclicGraphException(self)
		visited.add((self, requested_output))
		
		input, output, wildcards = self._expand_wildcards(requested_output)
		
		if output and (output, self) in jobs:
			return jobs[(output, self)]
		
		missing_input_exceptions = list()
		files_produced_with_error = set()
		todo = set()
		produced = dict()
		for rule, file in self._to_visit(input):
			try:
				job = rule.run(file, forceall = forceall, jobs = jobs, dryrun = dryrun, quiet = quiet, visited = set(visited))
				if file in produced:
					raise AmbiguousRuleException(produced[file], rule)
				if job.needrun:
					todo.add(job)
				produced[file] = rule
			except MissingInputException as ex:
				missing_input_exceptions.append(ex)
				files_produced_with_error.add(file)
		
		missing_input = self._get_missing_files(set(input) - produced.keys())
		if missing_input:
			raise MissingInputException(
				rule = self, 
				include = missing_input_exceptions, 
				files = set(missing_input) - files_produced_with_error
			)
		
		need_run = self._need_run(forcethis or forceall or todo, input, output)
		job = Job(
			self.workflow,
			rule = self, 
			message = self.get_message(input, output, wildcards),
			input = input,
			output = output,
			wildcards = Namedlist(fromdict = wildcards),
			depends = todo,
			dryrun = dryrun,
			needrun = need_run or quiet
		)
		jobs[(output, self)] = job
		return job

	def check(self):
		if self.output and not self.has_run():
			raise RuleException("Rule {} defines output but does not have a \"run\" definition.".format(self.name))

	def _need_run(self, force, input, output):
		""" Return True if rule needs to be run. """
		if self.has_run():
			if force:
				return True
			if self._has_missing_files(output):
				return True
			if not output:
				return True
			mintime = min(map(lambda f: os.stat(f).st_mtime, output))
			for i in input:
				if os.path.exists(i) and os.stat(i).st_mtime >= mintime: 
					
					return True
			return False
		return False

	def get_run(self):
		return self.workflow.get_run(self)

	def has_run(self):
		return self.workflow.has_run(self)

	def get_message(self, input, output, wildcards, showmessage = True):
		if self.message and showmessage:
			variables = dict(globals())
			variables.update(locals())
			return self.message.format(**variables)
		return "rule {}:\n\tinput: {}\n\toutput: {}\n".format(
			self.name, ", ".join(input), ", ".join(output))

	def is_producer(self, requested_output):
		"""
		Returns True if this rule is a producer of the requested output.
		"""
		for o in self.regex_output:
			match = re.match(o, requested_output)
			if match and len(match.group()) == len(requested_output):
				return True
		return False

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

	def __repr__(self):
		return self.name
