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
		self.parents = defaultdict(list)
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
		self.ignore_ambiguity = ignore_ambiguity
		
		targetjobs = map(partial(self.rule2job, self.rules), self.targetrules)
		self.update(targetjobs)
		self.check(children)
		
	def update(self, jobs):
		self.check(self.collect_children(jobs))
		
	def check(self, children, targetjobs):
		parents = self.parents
		if self.ignore_ambiguity:
			select_parent = itemgetter(0)
		else:
			def select_parent(jobs):
				jobs = sorted(jobs)
				if jobs[-1] > jobs[-2]:
					return job
				else:
					raise AmbiguousRuleException()
		
		# BFS
		queue = list(targetjobs)
		visited = set(queue)
		while queue:
			job = queue.pop(0)
			for file, jobs in children[job].items():
				job_ = select_parent(jobs)
				parents[job_] = job_
				if job_ not in visited:
					
	
	def collect_children(self, jobs):
		children = defaultdict(partial(defaultdict, list))
		file2jobs = partial(file2jobs, self.rules)
		
		# BFS
		queue = list(jobs) 
		visited = set()
		while queue:
			job = queue.pop(0)
			for file in job.input:
				for job_ in file2jobs(file):
					if job_ not in visited:
						queue.append(job_)
						visited.add(job_)
						children[job][file].append(job_)
						# TODO pumping lemma test!
		return children
	
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
