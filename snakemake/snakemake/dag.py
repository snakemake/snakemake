from collections import defaultdict
from itertools import chain, combinations
from functools import partial
from operator import itemgetter

from snakemake.jobs import Job

class DAG:
	def __init__(self, 
	             workflow,
	             targetfiles = None, 
	             targetrules = None,
	             forceall = False,
	             forcetargets = False, 
	             forcerules = None,
	             ignore_ambiguity = False):
		self.dependencies = defaultdict(partial(defaultdict, list))
		self.depending = defaultdict(partial(defaultdict, list))
		self.jobs = dict()
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
			self.select_dependency = self.select_dependency_ign_amb
		
		self.targetjobs = map(partial(self.rule2job, self.rules), self.targetrules)
		for job in self.targetjobs:
			self.update(job)
		
	def update(self, job, visited = None, parentmintime = None, skip_until_dynamic = False):
		if job in self.dependencies:
			return
		if visited is None:
			visited = set(job)
		dependencies = self.dependencies[job]
		potential_dependencies = self.collect_potential_dependencies(job).items()

		updated_input = list()		
		output_mintime = IOFile.mintime(job.output) or parentmintime
		skip_until_dynamic = not job.dynamic_output and (skip_until_dynamic or job.dynamic_input)
			
		for file, jobs in potential_dependencies.items():
			job_ = self.select_dependency(jobs, visited)
			self.update(job_, visited=set(visited), parentmintime=output_mintime, skip_until_dynamic=skip_until_dynamic)
			# TODO check for pumping up wildcards...
			if not skip_until_dynamic:
				dependencies[job_].append(file)
				self.depending[job_][job].append(file)
			if job_.needrun:
				updated_input.append(file)
				
		missing_input = job.rule.get_missing_files() - potential_dependencies.keys()
		if missing_input:
			raise MissingInputFileException()
		job.needrun, job.reason = self.needrun(self, job, updated_input, output_mintime)
		if needrun:
			self.update_needrun_temp(job)
	
	def needrun(self, job, updated_input, output_mintime):
		# forced?
		if job.rule in self.forcerules or (forcetargets and job in self.targetjobs):
			return True, "Forced execution"
			
		# missing output?
		output_missing = [f for f in job.output if not f.exists()]
		if output_missing:
			return True, "Missing output files: {}".format(", ".join(output_missing))
			
		# updated input?
		for f in job.input:
			if f.exists() and f.is_newer(output_mintime):
				updated_input.append(f)
		if updated_input:
			return True, "Updated input files: {}".format(", ".join(updated_input))
			
		return False, None
	
	def update_needrun_temp(self, job)
		for job_, files in self.dependencies[job].items():
			for f in files:
				if f in job.temp_output and not f.exists():
					job_.needrun, job.reason = True, "Missing temp files for rule {}".format(job)
					self.update_needrun_temp(job_)
					return
	
	def delete_job(self, job):
		for job_ in self.depending[job]:
			del self.dependencies[job_][job]
		del self.depending[job]
		for job_ in self.dependencies[job]:
			depending = self.depending[job_]
			del depending[job]
			if not depending:
				self.delete_job(job_)
		del self.dependencies[job]
	
	def collect_potential_dependencies(self, job):
		file2jobs = partial(file2jobs, self.rules)
		dependencies = defaultdict(list)
		for file in job.input:
			for job_ in file2jobs(file):
				dependencies[file].append(job_)
		return dependencies

	def select_dependency(self, jobs, visited):
		jobs = sorted(jobs)
		for i, job in reversed(enumerate(jobs)):
			if job not in visited and i > 0 and job > jobs[i-1]:
				return job
			else:
				raise AmbiguousRuleException()
		raise CyclicGraphException()
	
	def select_dependency_ign_amb(self, jobs, visited):
		for job in jobs:
			if job not in visited:
				return job
		raise CyclicGraphException()
	
	@lru_cache()
	@staticmethod
	def rule2job(rules, targetrule):
		return Job(rule=rule)
	
	@lru_cache()
	@staticmethod
	def file2jobs(rules, targetfile):
		for rule in rules:
			if rule.is_producer(file):
				yield Job(rule, targetfile=file)

