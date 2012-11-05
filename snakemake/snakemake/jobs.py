# -*- coding: utf-8 -*-

import os

from collections import defaultdict
from functools import lru_cache

from snakemake.utils import format

__author__ = "Johannes Köster"

class Job:
	def __init__(self, rule, targetfile = None):
		self.rule = rule
		self.targetfile = targetfile
		self.finished = False
		self.needrun = False
		self._hash = None
		
		self.input, self.output, self.log, self.wildcards = rule.expand_wildcards(self.targetfile)
		self.threads = rule.threads
		self.message = self._format_wildcards(rule.message) if rule.message else None
		self.shellcmd = self._format_wildcards(rule.shellcmd) if rule.shellcmd else None
		
		self.dynamic_output, self.dynamic_input, self.temp_output, self.protected_output = [set()]*4
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
		for i, f in enumerate(self.output):
			if f in self.dynamic_output:
				for f, _ in listfiles(self.rule.output[i]):
					yield IOFile(f)
			else:
				yield f
	
	@property
	def dynamic_wildcards(self):
		wildcards = defaultdict(set)
		for i, f in enumerate(self.output):
			if f in self.dynamic_output:
				for f, w in listfiles(self.rule.output[i]):
					for name, value in w.items():
						wildcards[name].add(value)
		return wildcards
	
	@property
	def missing_input(self):
		return set(f for f in self.input if not f.exists())
		
	def output_mintime(self):
		existing = [f.mtime() for f in self.output if f.exists()]
		if existing:
			return min(existing)
		return None
		
	def check_output(self):
		for f in self.expanded_output:
			if not f.exists():
				raise MissingOutputException("Output file {} not produced by rule {}.".format(f, self.rule.name), lineno = self.rule.lineno, snakefile = self.rule.snakefile)
	
	def protect_output(self):
		for f in self.expanded_output:
			if f in self.protected_output:
				f.protect()
	
	def cleanup(self):
		for f in self.output:
			if f.exists():
				f.remove()
	
	def _format_wildcards(self, string):
		return format(string, 
		              input=self.input, 
		              output=self.output, 
		              wildcards=self.wildcards, 
		              threads=self.threads, 
		              log=self.log)

	def __repr__(self):
		return self.rule.name
	
	def __eq__(self, other):
		return self.rule == other.rule and self.output == other.output
	
	def __hash__(self):
		if self._hash is None:
			self._hash = self.rule.__hash__()
			for o in self.output:
				self._hash ^= o.__hash__()
		return self._hash
