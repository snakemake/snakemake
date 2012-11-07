import textwrap
from collections import defaultdict
from itertools import chain, combinations, filterfalse
from functools import partial, lru_cache
from operator import itemgetter

from snakemake.io import IOFile
from snakemake.jobs import Job
from snakemake.exceptions import RuleException, MissingInputException, MissingRuleException, AmbiguousRuleException, CyclicGraphException

class DAG:
	def __init__(self, 
	             workflow,
	             targetfiles = None, 
	             targetrules = None,
	             forceall = False,
	             forcetargets = False, 
	             forcerules = None,
	             ignore_ambiguity = False):
		self.dependencies = defaultdict(partial(defaultdict, set))
		self.depending = defaultdict(partial(defaultdict, set))
		self._needrun = set()
		self._reason = dict()
		self._finished = set()
		self._dynamic = set()
		self._len = 0
		self.workflow = workflow
		self.rules = workflow.rules
		self.targetfiles = targetfiles
		self.targetrules = targetrules
		self.forcerules = set()
		self.forcefiles = set()
		if forceall:
			self.forcerules.update(self.rules)
		elif forcerules:
			self.forcerules.update(forcerules)
		if forcetargets:
			self.forcerules.update(targetrules)
			self.forcefiles.update(targetfiles)
		if ignore_ambiguity:
			self.select_dependency = self.select_dependency_ign_amb
		
		self.targetjobs = set(map(self.rule2job, self.targetrules))
		exceptions = defaultdict(list)		
		for file in self.targetfiles:
			try:
				for job in self.file2jobs(file):
					self.targetjobs.add(job)
			except MissingRuleException as ex:
				exceptions[file].append(ex)
		has_producer = set()
		
		for job in self.targetjobs:
			try:
				self.update(job)
				has_producer.add(job.targetfile)
			except MissingInputException as ex:
				exceptions[job.targetfile].append(ex)
		for file in has_producer:
			if file in exceptions:
				del exceptions[file]
		if exceptions:
			raise RuleException(include=chain(*exceptions.values()))
		
		for job in filter(lambda job: job.dynamic_output and not self.needrun(job), self.jobs):
			self.update_dynamic(job)
		self.update_needrun_temp()

	@property
	def jobs(self):
		for job in self.bfs(self.dependencies, *self.targetjobs):
			yield job
	
	@property
	def open_jobs(self):
		for job in filter(self.needrun, self.bfs(self.dependencies, *self.targetjobs, stop=self.finished)):
			yield job
			
	@property
	def ready_jobs(self):
		for job in filter(self.ready, self.open_jobs):
			yield job
			
	def ready(self, job):
		return all(map(lambda job: self.finished(job) or not self.needrun(job), self.dependencies[job])) and not self.dynamic(job) and not job.dynamic_input
	
	def needrun(self, job):
		return job in self._needrun
	
	def reason(self, job):
		return self._reason[job]
	
	def finished(self, job):
		return job in self._finished
	
	def dynamic(self, job):
		return job in self._dynamic
	
	def missing_temp(self, job):
		for job_, files in self.depending[job].items():
			if self.needrun(job_) and any(not f.exists() for f in files):
				return True
		return False
		
	def check_output(self, job):
		for f in job.expanded_output:
			if not f.exists():
				raise MissingOutputException("Output file {} not produced by rule {}.".format(f, job.rule.name), lineno = job.rule.lineno, snakefile = job.rule.snakefile)
	
	def handle_protected(self, job):
		for f in job.expanded_output:
			if f in job.protected_output:
				f.protect()
	
	def handle_temp(self, job):
		needed = lambda job, f: any(f in files for job_, files in self.depending[job].items() if not self.finished(job_))
		for job_ , files in self.dependencies[job].items():
			for f in files:
				if not needed(job_, f):
					f.remove()
		
	def update(self, job, visited = None, parentmintime = None, skip_until_dynamic = False, count = False, skipneedrun = False):
		if job in self.dependencies:
			return
		if visited is None:
			visited = set([job])
		dependencies = self.dependencies[job]
		potential_dependencies = self.collect_potential_dependencies(job).items()

		output_mintime = job.output_mintime or parentmintime
		#if job.rule.name == "annotate_clusters" and not skip_until_dynamic:
		#	import pdb; pdb.set_trace()
		skip_until_dynamic = skip_until_dynamic and not job.dynamic_output
		
		missing_input = job.missing_input
		exceptions = list()
		for file, jobs in potential_dependencies:
			try:
				job_ = self.select_dependency(file, jobs, visited)
				self.update(job_, visited=set(visited), parentmintime=output_mintime, skip_until_dynamic=skip_until_dynamic or file in job.dynamic_input)
				if file in missing_input:
					missing_input.remove(file)
				# TODO check for pumping up wildcards...
				dependencies[job_].add(file)
				self.depending[job_][job].add(file)
			except (RuleException, RuntimeError) as ex:
				if isinstance(ex, RuntimeError) and str(ex).startswith("maximum recursion depth exceeded"):
					raise RuleException("Maximum recursion depth exceeded. Maybe you have a cyclic dependency due to infinitely filled wildcards?\nProblematic input file:\n{}".format(file), lineno = job.rule.lineno, snakefile = job.rule.snakefile)
				if isinstance(ex, CyclicGraphException) and ex.file == file:
					raise ex
				exceptions.append(ex)
		if missing_input:
			raise MissingInputException(job.rule, missing_input)
		
		if skip_until_dynamic:
			self._dynamic.add(job)
		if not skipneedrun:
			needrun, reason = self._calc_needrun(job, output_mintime)
			if needrun:
				#print(job, job.output, reason, skip_until_dynamic)
				count = True
				self._needrun.add(job)
				self._reason[job] = reason
		if count:
			self._len += 1

	
	def _calc_needrun(self, job, output_mintime):
		#if job.rule.name == "list_clusters":
		#	import pdb; pdb.set_trace()
		# forced?
		if job.rule in self.forcerules or job.targetfile in self.forcefiles:
			return True, "Forced execution"
			
		# missing output?
		missing_output = job.missing_output
		if missing_output:
			return True, "Missing output files: {}".format(", ".join(missing_output))
		
		# updated input?
		updated_input = set(chain(*(f for job_, f in self.dependencies[job].items() if self.needrun(job_))))
		if output_mintime:
			for f in job.input:
				if f.exists() and f.is_newer(output_mintime):
					updated_input.add(f)
		if updated_input:
			return True, "Updated input files: {}".format(", ".join(updated_input))
			
		return False, None
	
	def finish(self, job, update_dynamic = True):
		self._finished.add(job)
		if update_dynamic:
			self.update_dynamic(job)
			self.update_needrun_temp()
	
	def update_dynamic(self, job):
		dynamic_wildcards = job.dynamic_wildcards
		if not dynamic_wildcards:
			# this happens e.g. in dryrun if output is not yet present
			return
		depending = list(filter(lambda job_: not self.finished(job_), self.bfs(self.depending, job)))
		job.rule.update_dynamic(dynamic_wildcards, input=False)
		newjob = Job(job.rule)
		self.replace_job(job, newjob)
		for job_ in depending:
			if job_.dynamic_input:
				if job_.rule.update_dynamic(dynamic_wildcards):
					if not self.dynamic(job_):
						#print(job_, job_.targetfile, job_.dynamic_output)
						newjob = Job(job_.rule, targetfile=job_.targetfile)
						self.replace_job(job_, newjob)
	
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
		if job in self._needrun:
			self._len -= 1
			self._needrun.remove(job)
			del self._reason[job]
		if job in self._finished:
			self._finished.remove(job)
		if job in self._dynamic:
			self._dynamic.remove(job)
	
	def replace_job(self, job, newjob):
		finished = self.finished(job)
		self.delete_job(job)
		self.update(newjob, count=finished, skipneedrun=finished)
		if job in self.targetjobs:
			self.targetjobs.remove(job)
			self.targetjobs.add(newjob)
	
	def update_needrun_temp(self):
		for job in filter(lambda job: job.temp_output, 
		                  self.bfs(self.dependencies, 
		                           *self.targetjobs, 
		                           stop=lambda job: self.finished(job) or not self.needrun(job) and not job.temp_output)):
			if self.missing_temp(job):
				self._needrun.add(job)
				self._reason[job] = "Missing temp file"
	
	def collect_potential_dependencies(self, job):
		dependencies = defaultdict(list)
		for file in job.input:
			try:
				for job_ in self.file2jobs(file):
					dependencies[file].append(job_)
			except MissingRuleException as ex:
				pass
		return dependencies

	def select_dependency(self, file, jobs, visited):
		jobs = sorted(jobs)
		for i, job in reversed(list(enumerate(jobs))):
			if job not in visited:
				if i == 0 or job > jobs[i-1]:
					return job
				else:
					raise AmbiguousRuleException(file, job.rule, jobs[i-1].rule)
		raise CyclicGraphException(job.rule, file)
	
	def select_dependency_ign_amb(self, file, jobs, visited):
		for job in jobs:
			if job not in visited:
				return job
		raise CyclicGraphException(job.rule, file)
	
	def bfs(self, direction, *jobs, stop=lambda job: False):
		queue = list(jobs)
		visited = set()
		while queue:
			job = queue.pop(0)
			if stop(job):
				# stop criterion reached for this node
				continue
			yield job
			for job_, _ in direction[job].items():
				if not job_ in visited:
					queue.append(job_)
					visited.add(job_)

	def new_wildcards(self, job):
		new_wildcards = set(job.wildcards.items())
		for job_ in self.dependencies[job]:
			if not new_wildcards:
				return set()
			for wildcard in job_.wildcards.items():
				new_wildcards.discard(wildcard)
		return new_wildcards

	def __str__(self):
		jobid = dict((job, i) for i, job in enumerate(self.jobs))
		nodes, edges = list(), list()
		for job in self.jobs:
			label = "\\n".join([job.rule.name] + list(map(": ".join, self.new_wildcards(job))))
			style = ""
			if not self.needrun(job):
				style = ',style="rounded,dashed"'
			if self.dynamic(job) or job.dynamic_input:
				style = ',style="rounded,dotted"'
			nodes.append('\t{}[label = "{}"{}];'.format(jobid[job], label, style))
			for job_ in self.dependencies[job]:
				edges.append("\t{} -> {};".format(jobid[job_], jobid[job]))
		return textwrap.dedent("""\
		                    digraph snakemake_dag {{
		                    	node[shape=box,style=rounded];
		                    {nodes}
		                    {edges}
		                    }}\
		                    """).format(nodes="\n".join(nodes), edges="\n".join(edges))

	def __len__(self):
		return self._len
	
	def rule2job(self, targetrule):
		return Job(rule=targetrule)
	
	def file2jobs(self, targetfile):
		jobs = list()
		for rule in self.rules:
			if rule.is_producer(targetfile):
				jobs.append(Job(rule, targetfile=targetfile))
		if not jobs:
			raise MissingRuleException(targetfile)
		return jobs
