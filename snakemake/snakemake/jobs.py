# -*- coding: utf-8 -*-

from snakemake.utils import format

__author__ = "Johannes Köster"

class Job:
	def __init__(self, rule, targetfile = None):
		self.rule = rule
		self.targetfile = targetfile
		self.finished = False
		
		self.input, self.output, self.log, self.wildcards = rule.expand_wildcards(self.targetfile)
		self.message = self._format_wildcards(rule.message)
		self.shellcmd = self._format_wildcards(rule.shellcmd) if rule.shellcmd else None
		
		self.dynamic_output, self.temp_output, self.protected_output = set(), set(), set()
		for i, f in self.output:
			f_ = self.rule.output[i]
			if f_ in self.rule.dynamic_output:
				self.dynamic_output.add(f)
			if f_ in self.rule.temp_output:
				self.temp_output.add(f)
			if f_ in self.rule.protected_output:
				self.protected_output.add(f)
	
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
	
	def _format_wildcards(self, string):
		return format(string, 
		              input=self.input, 
		              output=self.output, 
		              wildcards=self.wildcards, 
		              threads=self.threads, 
		              log=self.log)

	def __str__(self):
		return self.rule.name
	
	@lru_cache()			
	def __eq__(self, other):
		return self.rule == other.rule and self.output == other.output
	
	@lru_cache()
	def __hash__(self):
		h = self.rule.__hash__()
		for o in self.output:
			h ^= o.__hash__()
		return h
