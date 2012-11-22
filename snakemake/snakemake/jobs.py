# -*- coding: utf-8 -*-

import os, sys

from collections import defaultdict
from itertools import chain
from functools import lru_cache

from snakemake.io import IOFile
from snakemake.utils import format, listfiles
from snakemake.exceptions import MissingOutputException, RuleException

__author__ = "Johannes Köster"

class Job:
	HIGHEST_PRIORITY = sys.maxsize

	def __init__(self, rule, targetfile = None):
		self.rule = rule
		self.targetfile = targetfile
		self._hash = None
		
		self.input, self.output, self.log, self.wildcards = rule.expand_wildcards(self.targetfile)
		self.threads = rule.threads
		self.priority = rule.priority
		try:
			self.message = self._format_wildcards(rule.message) if rule.message else None
			self.shellcmd = self._format_wildcards(rule.shellcmd) if rule.shellcmd else None
		except AttributeError as ex:
			raise RuleException(str(ex), rule=self.rule)
		except KeyError as ex:
			raise RuleException("Unknown variable in message of shell command: {}".format(str(ex)), rule=self.rule)
		
		
		self.dynamic_output, self.dynamic_input, self.temp_output, self.protected_output = set(), set(), set(), set()
		for f, f_ in zip(self.output, self.rule.output):
			if f_ in self.rule.dynamic_output:
				self.dynamic_output.add(f)
			if f_ in self.rule.temp_output:
				self.temp_output.add(f)
			if f_ in self.rule.protected_output:
				self.protected_output.add(f)
		for f, f_ in zip(self.input, self.rule.input):
			if f_ in self.rule.dynamic_input:
				self.dynamic_input.add(f)
	
	@property
	def expanded_output(self):
		for f, f_ in zip(self.output, self.rule.output):
			if f in self.dynamic_output:
				expansion = self.expand_dynamic(f_)
				if not expansion:
					yield f_
				for f, _ in expansion:
					yield IOFile(f, self.rule)
			else:
				yield f
	
	@property
	def dynamic_wildcards(self):
		wildcards = defaultdict(set)
		for f, f_ in zip(self.output, self.rule.output):
			if f in self.dynamic_output:
				for f, w in self.expand_dynamic(f_):
					for name, value in w.items():
						wildcards[name].add(value)
		return wildcards
	
	@property
	def missing_input(self):
		return set(f for f in self.input if not f.exists)
	
	@property
	def output_mintime(self):
		existing = [f.mtime for f in self.expanded_output if f.exists]
		if existing:
			return min(existing)
		return None
	
	def missing_output(self, requested = None):
		if requested is None:
			requested = set(self.output)
		files = set()
		
		for f, f_ in zip(self.output, self.rule.output):
			if f in requested:
				if f in self.dynamic_output:
					if not self.expand_dynamic(f_):
						files.add("{} (dynamic)".format(f_))
				elif not f.exists:
					files.add(f)
		return files
	
	def prepare(self):
		protected = list(filter(lambda f: f.protected, self.expanded_output))
		if protected:
			raise ProtectedOutputException(self.rule, protected)
		if self.dynamic_output:
			for f, _ in chain(*map(self.expand_dynamic, self.rule.dynamic_output)):
				os.remove(f)
		for f, f_ in zip(self.output, self.rule.output):
			if f in self.dynamic_output:
				f_.prepare()
			else:
				f.prepare()
		if self.log:
			self.log.prepare()
	
	def cleanup(self):
		for f in self.output:
			if f.exists:
				f.remove()
	
	def _format_wildcards(self, string):
		return format(string, 
		              input=self.input, 
		              output=self.output, 
		              wildcards=self.wildcards, 
		              threads=self.threads, 
		              log=self.log, **self.rule.workflow.globals)

	def __repr__(self):
		return self.rule.name
	
	def __eq__(self, other):
		if other is None:
			return False
		return self.rule == other.rule and self.output == other.output

	def __lt__(self, other):
		return self.rule.__lt__(other.rule)

	def __gt__(self, other):
		return self.rule.__gt__(other.rule)
	
	def __hash__(self):
		if self._hash is None:
			self._hash = self.rule.__hash__()
			for o in self.output:
				self._hash ^= o.__hash__()
		return self._hash
	
	@staticmethod
	def expand_dynamic(pattern):
		return list(listfiles(pattern))

class Reason:
	def __init__(self):
		self.updated_input = set()
		self.updated_input_run = set()
		self.missing_output = set()
		self.forced = False

	def __str__(self):
		if self.forced:
			return "Forced execution"
		if self.missing_output:
			return "Missing output files: {}".format(", ".join(self.missing_output))
		if self.updated_input:
			return "Updated input files: {}".format(", ".join(self.updated_input))
		if self.updated_input_run:
			return "This run will update input files: {}".format(", ".join(self.updated_input_run))
		return ""

	def __bool__(self):
		return bool(self.updated_input or self.missing_output or self.forced or self.updated_input_run)
