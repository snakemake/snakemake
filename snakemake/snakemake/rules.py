# -*- coding: utf-8 -*-

import os, re, sys, inspect, sre_constants
from collections import defaultdict

from snakemake.io import IOFile, protected, temp, dynamic, Namedlist, expand
from snakemake.exceptions import RuleException

__author__ = "Johannes Köster"

class Rule:
	def __init__(self, name, workflow, lineno = None, snakefile = None):
		"""
		Create a rule
		
		Arguments
		name -- the name of the rule
		"""
		self.name = name
		self.workflow = workflow
		self.docstring = None
		self.message = None
		self._input = Namedlist()
		self._output = Namedlist()
		self.dynamic_output = set()
		self.dynamic_input = set()
		self.temp_output = set()
		self.protected_output = set()
		self.threads = 1
		self.priority = 1
		self._log = None
		self.wildcard_names = set()
		self.lineno = lineno
		self.snakefile = snakefile
		self.run_func = None
		self.shellcmd = None

	def update_dynamic(self, wildcards, input=True):
		io, dynamic_io = (self.input, self.dynamic_input) if input else (self.output, self.dynamic_output)
		expansion = defaultdict(list)
		for i, f in enumerate(io):
			if f in dynamic_io:
				try:
					for e in reversed(expand(f, zip, **wildcards)):
						expansion[i].append(IOFile(e, rule=self))
				except KeyError:
					return False
		# replace the dynamic input with the expanded files
		for i, e in reversed(list(expansion.items())):
			dynamic_io.remove(io[i])
			io.insert_items(i, e)
		if not input:
			self.wildcard_names.clear()
		return True
		

	def has_wildcards(self):
		"""
		Return True if rule contains wildcards.
		"""
		return bool(self.wildcard_names)
	
	@property
	def log(self):
		return self._log
	
	@log.setter
	def log(self, log):
		self._log = IOFile(log, rule=self)

	@property
	def input(self):
		return self._input

	def set_input(self, *input, **kwinput):
		"""
		Add a list of input files. Recursive lists are flattened.
		
		Arguments
		input -- the list of input files
		"""
		for item in input:
			self._set_inoutput_item(item)
		for name, item in kwinput.items():
			self._set_inoutput_item(item, name = name)

	@property
	def output(self):
		return self._output

	def set_output(self, *output, **kwoutput):
		"""
		Add a list of output files. Recursive lists are flattened.
		
		Arguments
		output -- the list of output files
		"""
		for item in output:
			self._set_inoutput_item(item, output = True)
		for name, item in kwoutput.items():
			self._set_inoutput_item(item, output = True, name = name)
		
		for item in self.output:
			wildcards = item.get_wildcard_names()
			if self.wildcard_names:
				if self.wildcard_names != wildcards:
					raise SyntaxError("Not all output files of rule {} contain the same wildcards. ".format(self.name))
			else:
				self.wildcard_names = wildcards
	
	def _set_inoutput_item(self, item, output = False, name=None):
		"""
		Set an item to be input or output.
		
		Arguments
		item     -- the item
		inoutput -- either a Namedlist of input or output items
		name     -- an optional name for the item
		"""
		inoutput = self.output if output else self.input
		if inspect.isfunction(item) and output:
			raise SyntaxError("Only input files can be specified as functions")
		if isinstance(item, str) or type(item).__name__ == "function":
			_item = IOFile(item, rule=self)
			if isinstance(item, temp):
				if not output:
					raise SyntaxError("Only output files may be temporary")
				self.temp_output.add(_item)
			if isinstance(item, protected):
				if not output:
					raise SyntaxError("Only output files may be protected")
				self.protected_output.add(_item)
			if isinstance(item, dynamic):
				if len(_item.get_wildcard_names()) > 1:
					raise SyntaxError("Dynamic files may not contain more than one wildcard.")
				if output:
					self.dynamic_output.add(_item)
				else:
					self.dynamic_input.add(_item)
			inoutput.append(_item)
			if name:
				inoutput.add_name(name)
		else:
			try:
				for i in item:
					self._set_inoutput_item(i, output = output)
			except TypeError:
				raise SyntaxError("Input and output files must be specified as strings.")
		
	def expand_wildcards(self, requested_output):
		""" Expand wildcards depending on the requested output. """
		wildcards = dict()
		if requested_output:
			wildcards = self.get_wildcards(requested_output)
			missing_wildcards = set(wildcards.keys()) - self.wildcard_names 
		else:
			missing_wildcards = self.wildcard_names
		
		if missing_wildcards:
			raise RuleException("Could not resolve wildcards in rule {}:\n{}".format(self.name, "\n".join(self.wildcard_names)), lineno = self.lineno, snakefile = self.snakefile)

		try:
			input = Namedlist()
			for f in self.input:
				if f in self.dynamic_input:
					input.append(f.fill_wildcards())
				else:
					input.append(f.apply_wildcards(wildcards))
			output = Namedlist(o.apply_wildcards(wildcards) for o in self.output)
			input.take_names(self.input.get_names())
			output.take_names(self.output.get_names())
			log = self.log.apply_wildcards(wildcards) if self.log else None
			return input, output, log, Namedlist(fromdict = wildcards)
		except KeyError as ex:
			# this can only happen if an input file contains an unresolved wildcard.
			raise RuleException("Wildcards in input or log file of rule {} do not appear in output files:\n{}".format(self, str(ex)), lineno = self.lineno, snakefile = self.snakefile)

	def is_producer(self, requested_output):
		"""
		Returns True if this rule is a producer of the requested output.
		"""
		try:
			for o in self.output:
				match = o.match(requested_output)
				if match and len(match.group()) == len(requested_output):
					return True
			return False
		except sre_constants.error as ex:
			raise IOFileException("{} in wildcard statement".format(ex), snakefile=self.snakefile, lineno=self.lineno)

	def get_wildcards(self, requested_output):
		"""
		Update the given wildcard dictionary by matching regular expression output files to the requested concrete ones.
		
		Arguments
		wildcards -- a dictionary of wildcards
		requested_output -- a concrete filepath
		"""
		bestmatchlen = 0
		bestmatch = None
		bestmatch_output = None
		for i, o in enumerate(self.output):
			match = o.match(requested_output)
			if match:
				l = self.get_wildcard_len(match.groupdict())
				if not bestmatch or bestmatchlen > l:
					bestmatch = match.groupdict()
					bestmatchlen = l
					bestmatch_output = self.output[i]
		return bestmatch
	
	@staticmethod
	def get_wildcard_len(wildcards):
		"""
		Return the length of the given wildcard values.
		
		Arguments
		wildcards -- a dict of wildcards
		"""
		return sum(map(len, wildcards.values()))

	def __lt__(self, rule):
		comp = self.workflow.ruleorder.compare(self.name, rule.name)
		return comp < 0

	def __gt__(self, rule):
		comp = self.workflow.ruleorder.compare(self.name, rule.name)
		return comp > 0

	def __str__(self):
		return self.name


class Ruleorder:
	def __init__(self):
		self.order = list()

	def add(self, *rulenames):
		"""
		Records the order of given rules as rule1 > rule2 > rule3, ...
		"""
		self.order.append(list(rulenames))

	def compare(self, rule1name, rule2name):
		"""
		Return whether rule2 has a higher priority that rule1.
		"""
		# try the last clause first, i.e. clauses added later overwrite those before.
		for clause in reversed(self.order):
			try:
				i = clause.index(rule1name)
				j = clause.index(rule2name)
				# rules with higher priority should have a smaller index
				comp = j - i
				if comp < 0: 
					comp = -1
				elif comp > 0:
					comp = 1
				return comp
			except ValueError:
				pass
		return 0
