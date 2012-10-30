# -*- coding: utf-8 -*-

import signal
import sys, time, os, threading, multiprocessing
from itertools import chain, filterfalse
from collections import defaultdict
from snakemake.exceptions import TerminatedException, MissingOutputException, RuleException, \
	ClusterJobException, print_exception, format_error, get_exception_origin
from snakemake.shell import shell
from snakemake.io import IOFile, temp, protected, expand, touch, remove
from snakemake.utils import listfiles
from snakemake.logging import logger
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

if os.name == "posix":
	from multiprocessing import Event
	PoolExecutor = ProcessPoolExecutor
else:
	from threading import Event
	PoolExecutor = ThreadPoolExecutor

__author__ = "Johannes Köster"

def run_wrapper(run, rulename, ruledesc, input, output, wildcards, 
		threads, log, rowmaps, rulelineno, rulesnakefile):
	"""
	Wrapper around the run method that handles directory creation and
	output file deletion on error.
	
	Arguments
	run -- the run method
	input -- list of input files
	output -- list of output files
	wildcards -- so far processed wildcards
	"""
	logger.info(ruledesc)

			
	t0 = time.time()
	try:
		# execute the actual run method.
		run(input, output, wildcards, threads, log)
		# finish all spawned shells.
		shell.join_all()
		t1 = time.time()
		return t0, t1
	except (Exception, BaseException) as ex:
		# this ensures that exception can be re-raised in the parent thread
		lineno, file = get_exception_origin(ex, rowmaps)
		raise RuleException(format_error(ex, lineno, rowmaps=rowmaps, snakefile=file), )



class Job:
	def __init__(self, rule, targetfile = None):
		self.rule = rule
		self.targetfile = targetfile
	
	def init(self):
		self.input, self.output, self.log, self.wildcards = rule.expand_wildcards(self.targetfile)
	
	def needrun(self):
		pass


def print_job_dag(jobs):
	print("digraph snakemake_dag {")
	for job in jobs:
		for edge in job.dot():
			print("\t" + edge)
	print("}")



class Job:
	count = 0

	@staticmethod
	def cleanup_unfinished(jobs):
		for job in jobs:
			job.cleanup()

	def __init__(self, workflow, rule = None, message = None, reason = None,
			input = None, output = None, wildcards = None, shellcmd = None,
			threads = 1, log = None, depends = set(), dryrun = False, quiet = False,
			touch = False, needrun = True, pseudo = False, visited = None, dynamic_output = False, forced = False, forcerules = None, forceall = False):
		self.workflow = workflow
		self.scheduler = None
		self.rule = rule
		self.message = message
		self.reason = reason
		self.input = input
		self.output = output
		self.wildcards = wildcards
		self.threads = threads
		self.log = log
		self.dryrun = dryrun
		self.quiet = quiet
		self.touch = touch
		self.shellcmd = shellcmd
		self.needrun = needrun
		self.pseudo = pseudo
		self.forced = forced
		self.forcerules = forcerules
		self.forceall = forceall
		self.dynamic_output = dynamic_output
		self.depends = set(depends)
		self.depending = set()
		self.is_finished = False
		self._callbacks = list()
		self._error_callbacks = list()
		self.jobid = Job.count
		self.visited = visited
		self.ignore = False
		Job.count += 1
		for other in self.depends:
			other.depending.add(self)

	def all_jobs(self):
		yield self
		for job in self.descendants():
			yield job

	def descendants(self):
		for job in self.depends:
			yield job
			for j in job.descendants():
				yield j

	def ancestors(self):
		queue = [job for job in self.depending]
		visited = set(queue)
		while queue:
			job = queue.pop(0)
			for j in job.depending:
				if j not in visited:
					queue.append(j)
					visited.add(j)
			yield job
	
	def get_message(self):
		msg = ""
		if not self.quiet:
			if self.message:
				msg += self.message.format(input=self.input, output=self.output, wildcards=self.wildcards, threads=self.threads, log=self.log, **globals())
			else:
				def showtype(orig_iofiles, iofiles):
					for i, iofile in enumerate(iofiles):
						if self.rule.is_dynamic(iofile):
							f = str(orig_iofiles[i])
							f += " (dynamic)"
						else:
							f = str(iofile)
						if iofile.is_temp():
							f += " (temporary)"
						if iofile.is_protected():
							f += " (protected)"
						yield f
				
				msg += "rule " + self.rule.name
				if self.input or self.output:
					msg += ":"
				if self.input:
					msg += "\n\tinput: {}".format(", ".join(showtype(self.rule.input, self.input)))
				if self.output:
					msg += "\n\toutput: {}".format(", ".join(showtype(self.rule.output, self.output)))
				if self.reason:
					msg += "\n\t{}".format(self.reason)
		if self.shellcmd:
			if not self.quiet:
				msg += "\n"
			msg += self.shellcmd.format(input=self.input, output=self.output, wildcards=self.wildcards, threads=self.threads, log=self.log, **globals())
		return msg
		
	def print_message(self):
		logger.info(self.get_message())
		
	def run(self, run_func):
		if not self.needrun or self.pseudo or self.ignore:
			self.finished()
		elif self.dryrun:
			self.print_message()
			self.finished()
		elif self.touch:
			logger.info(self.message)
			for o in self.output:
				if self.rule.is_dynamic(o):
					for f, _ in listfiles(o):
						touch(f)
				else:
					o.touch(self.rule.name, self.rule.lineno, self.rule.snakefile)
			# sleep shortly to ensure that output files of different rules 
			# are not touched at the same time.
			time.sleep(0.1)
			self.finished()
		else:
			for o in self.output:
				if self.rule.is_dynamic(o):
					for f, _ in listfiles(o):
						try:
							IOFile(f).remove()
						except OSError:
							raise RuleException("Could not remove dynamic output file {}.".format(f), lineno=self.rule.lineno, snakefile=self.rule.snakefile)
				# TODO what if a directory inside o is dynamic?
				o.prepare()
			if self.log:
				self.log.prepare()

			run_func(self)
	
	def get_run_args(self):
		return (self.rule.get_run(), self.rule.name, self.get_message(), 
			self.input, self.output, self.wildcards, self.threads, self.log, 
			self.workflow.rowmaps, self.rule.lineno, self.rule.snakefile)
	
	def add_callback(self, callback):
		""" Add a callback that is invoked when job is finished. """
		self._callbacks.append(callback)

	def add_error_callback(self, callback):
		self._error_callbacks.append(callback)
	
	def finished(self, future = None):
		""" Set job to be finished. """
		self.is_finished = True
		if not self.ignore:
			if self.needrun and not self.pseudo:
				try:
					if future:
						ex = future.exception()
						if ex:
							raise ex

					if not self.dryrun:
						# check the produced files
						for o in self.output:
							if not self.rule.is_dynamic(o):
								o.created(self.rule.name, self.rule.lineno, self.rule.snakefile)
						for f in self.input:
							f.used()
				except (Exception, BaseException) as ex:
					# in case of an error, execute all callbacks and delete output
					print_exception(ex, self.workflow.rowmaps)
					self.cleanup()
					for callback in self._error_callbacks:
						callback()
					return

				if not self.dryrun:
					self.workflow.jobcounter.done()
					if not self.quiet:
						logger.info(self.workflow.jobcounter)
					if not future is None:
						self.workflow.report_runtime(self.rule, future.result())

			for other in self.depending:
				other.depends.remove(self)

			if not self.dryrun and self.dynamic_output:
				self.handle_dynamic_output()

		for callback in self._callbacks:
			callback(self)

	def cleanup(self):
		if not self.is_finished:
			for o in self.output:
				if self.rule.is_dynamic(o):
					for f, _ in listfiles(o):
						remove(f)
				else:
					o.remove()
				
	def handle_dynamic_output(self):
		wildcard_expansion = defaultdict(set)
		for i, o in enumerate(self.output):
			if self.rule.is_dynamic(o):
				for f, wildcards in listfiles(self.rule.output[i]):
					for name, value in wildcards.items():
						wildcard_expansion[name].add(value)
		#import pdb; pdb.set_trace()
		# determine jobs to add
		new_jobs = set()
		dynamic = 0
		jobs = {(self.output, self.rule): self}
		for job in self.scheduler.get_jobs():
			if not job.pseudo:
				jobs[(job.output, job.rule)] = job
		for job in self.ancestors():
			j = job.handle_dynamic_input(wildcard_expansion, jobs)
			if j:
				new_jobs.update(j.all_jobs())
				dynamic += 1
		# remove this job from the DAG as it would induce the dynamic loop again
		if self in new_jobs:
			new_jobs.remove(self)
			for job in list(self.depending):
				if self in job.depends:
					job.depends.remove(self)
					self.depending.remove(job)

		# calculate how many jobs have to be added
		n = len(new_jobs) - dynamic
		if n:
			logger.warning("Dynamically adding {} new jobs".format(n))
		self.scheduler.add_jobs(new_jobs)
		self.workflow.jobcounter.count += n

	def handle_dynamic_input(self, wildcard_expansion, jobs):
		expansion = defaultdict(list)
		for i, f in enumerate(self.rule.input):
			if self.rule.is_dynamic(f):
				try:
					for e in reversed(expand(f, zip, **wildcard_expansion)):
						expansion[i].append(IOFile.create(e, temp = f.is_temp(), protected = f.is_protected()))
				except Exception as ex:
					# keep the file if expansion fails
					return
		if not expansion:
			return
		# replace the dynamic input with the expanded files
		for i, e in reversed(list(expansion.items())):
			self.rule.set_dynamic(self.rule.input[i], False)
			#self.rule.input[i:i+1] = e
			self.rule.input.insert_items(i, e)
		if self.pseudo:
			return
		try:
			try:
				del jobs[self.output, self.rule]
			except KeyError:
				# job is not registered yet
				pass
			
			non_dynamic_output = [o for o in self.output if not self.rule.is_dynamic(o)]
			job = self.rule.run(non_dynamic_output[0] if non_dynamic_output else None, dryrun = self.dryrun, jobs=jobs, forcethis=self.forced, give_reason=self.reason, forceall=self.forceall, forcerules=self.forcerules)

			# remove current job from DAG
			for j in self.depends:
				j.depending.remove(self)
			for j in list(self.depending):
				j.depends.remove(self)
				j.depends.add(job)
				job.depending.add(j)
			self.depends = list()
			self.ignore = True

			return job
		except RuleException as ex:
			# there seem to be missing files, so ignore this
			pass
		

	def dot(self):
		label = self.rule.name
		new_wildcards = self.new_wildcards()
		if not self.depends or new_wildcards:
			for wildcard, value in new_wildcards:
				label += "\\n{}: {}".format(wildcard, value)
		edges = ("{} -> {};".format(j.jobid, self.jobid) 
			for j in self.depends if j.needrun)
		node = ('{}[label = "{}"];'.format(self.jobid, label),)
		return chain(node, edges)

	def new_wildcards(self):
		new_wildcards = set(self.wildcards.items())
		for job in self.depends:
			if not new_wildcards:
				return set()
			for wildcard in job.wildcards.items():
				new_wildcards.discard(wildcard)
		return new_wildcards

	def __repr__(self):
		return self.rule.name

