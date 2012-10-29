from collections import defaultdict
from itertools import chain, combinations
from functools import partial
from operator import itemgetter


class DAG:
	def __init__(self, 
	             workflow,
	             targetfiles = None, 
	             targetrules = None,
	             forceall = False,
	             forcetargets = False, 
	             forcerules = None,
	             ignore_ambiguity = False):
		self.dependencies = defaultdict(list)
		self.workflow = workflow
		self.rules = workflow.rules
		self.targetfiles = targetfiles
		self.targetrules = targetrules
		self.forcerules = set()
		if forceall:
			self.forcerules = self.rules
		elif forcerules:
			self.forcerules = forcerules
		if forcetargets:
			self.forcerules.update(targetrules)
		if ignore_ambiguity:
			self.select_parent = self.select_parent_ign_amb
		
		targetjobs = map(partial(self.rule2job, self.rules), self.targetrules)
		for job in targetjobs:
			self.update(job)
		
	def update(self, job, visited = None):
		if visited is None:
			visited = set(job)
		dependencies = self.dependencies[job]
		potential_dependencies = self.collect_dependencies(job).items()
		for file, jobs in potential_dependencies:
			job_ = self.select_dependencies(jobs, visited)
			self.update(job_, visited=set(visited))
			dependencies.append(job_)
		missing_input = #TODO missing input - potential_dependencies.keys()
		if missing_input:
			raise MissingInputFileException()
		
		
	
	def collect_dependencies(self, job):
		file2jobs = partial(file2jobs, self.rules)
		dependencies = defaultdict(list)
		for file in job.input:
			for job_ in file2jobs(file):
				dependencies[file].append(job_)
		return dependencies
	
	@lru_cache()
	@staticmethod
	def rule2job(rules, targetrule):
		return Job(rule=rule, *rule.expand_wildcards(None))
	
	@lru_cache()
	@staticmethod
	def file2jobs(rules, targetfile):
		for rule in rules:
			if rule.is_producer(file):
				yield Job(rule, *rule.expand_wildcards(file))
	
	def select_parent(self, jobs, visited):
		jobs = sorted(jobs)
		for i, job in reversed(enumerate(jobs)):
			if job not in visited and i > 0 and job > jobs[i-1]:
				return job
			else:
				raise AmbiguousRuleException()
		raise CyclicGraphException()
	
	def select_parent_ign_amb(self, jobs, visited):
		for job in jobs:
			if job not in visited:
				return job
		raise CyclicGraphException()
