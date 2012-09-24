# -*- coding: utf-8 -*-

import os, re, sys
import sre_constants
from collections import defaultdict
from snakemake.logging import logger
from snakemake.jobs import Job
from snakemake.io import IOFile, protected, temp, dynamic, Namedlist
from snakemake.utils import listfiles
from snakemake.exceptions import MissingInputException, AmbiguousRuleException, CyclicGraphException, RuleException, ProtectedOutputException, IOFileException

__author__ = "Johannes Köster"

class Rule:
	def __init__(self, *args, lineno = None, snakefile = None):
		"""
		Create a rule
		
		Arguments
		name -- the name of the rule
		"""
		if len(args) == 2:
			name, workflow = args
			self.name = name
			self.message = None
			self.input = Namedlist()
			self.output = Namedlist()
			self.dynamic = dict()
			self.threads = 1
			self.wildcard_names = set()
			self.workflow = workflow
			self.lineno = lineno
			self.snakefile = snakefile
			self.run_func = None
		elif len(args) == 1:
			other = args[0]
			self.name = other.name
			self.message = other.message
			self.input = other.input
			self.output = other.output
			self.dynamic = other.dynamic
			self.threads = other.threads
			self.wildcard_names = other.wildcard_names
			self.workflow = other.workflow
			self.lineno = other.lineno
			self.snakefile = other.snakefile
			self.run_func = other.run_func
		else:
			raise ValueError("Rule expects either 1 or 2 positional args.")


	def has_wildcards(self):
		"""
		Return True if rule contains wildcards.
		"""
		return bool(self.wildcard_names)

	def is_dynamic(self, file):
		return file in self.dynamic

	def set_dynamic(self, file, dynamic):
		if dynamic:
			# create a bipartite graph that assigns a concrete to each dynamic file
			concrete = file.fill_wildcards()
			self.dynamic[file] = concrete
			self.dynamic[concrete] = file
		else:
			concrete = self.dynamic[file]
			del self.dynamic[file]
			del self.dynamic[concrete]
	
	def get_non_dynamic_output(self, output):
		return [f for i, f in enumerate(output) if not self.is_dynamic(self.output[i])]	

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
		if type(item).__name__ == "function" and output:
			raise SyntaxError("Only input files can be specified as functions")
		try:
			_item = IOFile.create(item, temp = isinstance(item, temp), protected = isinstance(item, protected))
			inoutput.append(_item)
			if name:
				inoutput.add_name(name)
			if isinstance(item, dynamic):
				if len(_item.get_wildcard_names()) > 1:
					raise SyntaxError("Dynamic files may not contain more than one wildcard.")
				self.set_dynamic(_item, True)
		except ValueError:
			try:
				for i in item:
					self._set_inoutput_item(i, output = output)
			except TypeError:
				raise SyntaxError("Input and output files must be specified as strings.")

	def set_message(self, message):
		"""
		Set the message that is displayed when rule is executed.
		"""
		self.message = message
	
	def set_threads(self, threads):
		self.threads = threads
		
	def _expand_wildcards(self, requested_output):
		""" Expand wildcards depending on the requested output. """
		wildcards = dict()
		if requested_output:
			wildcards, matching_output = self.get_wildcards(requested_output)
			missing_wildcards = set(wildcards.keys()) - self.wildcard_names 
		else:
			missing_wildcards = self.wildcard_names
			matching_output = None
		
		if missing_wildcards:
			raise RuleException("Could not resolve wildcards in rule {}:\n{}".format(self.name, "\n".join(self.wildcard_names)), lineno = self.lineno, snakefile = self.snakefile)

		try:
			input = Namedlist()
			for i in self.input:
				if self.is_dynamic(i):
					input.append(self.dynamic[i])
				else:
					input.append(i.apply_wildcards(wildcards))
			output = Namedlist(o.apply_wildcards(wildcards) for o in self.output)
			input.take_names(self.input.get_names())
			output.take_names(self.output.get_names())
			return input, output, wildcards, matching_output
		except KeyError as ex:
			# this can only happen if an input file contains an unresolved wildcard.
			raise RuleException("Wildcards in input file of rule {} do not appear in output files:\n{}".format(self, str(ex)), lineno = self.lineno, snakefile = self.snakefile)

	@staticmethod
	def _get_missing_files(files):
		""" Return the tuple of files that are missing form the given ones. """
		return tuple(f for f in files if not os.path.exists(f))
	
	@staticmethod
	def _has_missing_files(files, requested):
		""" Return True if any of the given files does not exist. """
		for f in files:
			if (requested == None or f.get_file() in requested) and not os.path.exists(f) and not f.is_temp():
				return True
		return False
	
	def run(self, requested_output = None, forceall = False, forcethis = False, 
	        give_reason = False, jobs = None, dryrun = False, touch = False, 
	        quiet = False, visited = None, parentmintime = None, 
	        ignore_ambiguity = False, skip_until_dynamic = False):
		"""
		Run the rule.
		
		Arguments
		requested_output -- the optional requested output file
		forceall         -- whether all required rules shall be executed, even if the files already exist
		forcethis        -- whether this rule shall be executed, even if the files already exist
		jobs             -- dictionary containing all jobs currently queued
		dryrun           -- whether rule execution shall be only simulated
		visited          -- set of already visited pairs of rules and requested output
		"""
		if jobs is None:
			jobs = dict()
		if visited is None:
			visited = set()

		if (self, requested_output) in visited:
			raise CyclicGraphException(self, lineno = self.lineno, snakefile = self.snakefile)
		visited.add((self, requested_output))
		
		input, output, wildcards, matching_output = self._expand_wildcards(requested_output)
		
		skip_until_dynamic = skip_until_dynamic and not self.is_dynamic(matching_output)
		pseudo = skip_until_dynamic
	
		output_mintime = IOFile.mintime(output) or parentmintime
		
		exceptions = defaultdict(list)
		todo = set()
		produced = dict()

		if (output, self) in jobs:
			job = jobs[(output, self)]
			if not job.needrun:
				# update the job if it needs to run due to the new requested file
				needrun, reason = self._need_run(False, todo, input, output, parentmintime, requested_output, pseudo)
				job.needrun = needrun
				job.message = self.get_message(input, output, wildcards,
				                               reason if give_reason else None)
			return job
		
		for rule, file in self.workflow.get_producers(input, exclude=self):
			try:
				job = rule.run(
					file, 
					forceall = forceall, 
					jobs = jobs, 
					dryrun = dryrun, 
					give_reason = give_reason,
					touch = touch, 
					quiet = quiet, 
					visited = set(visited), 
					parentmintime = output_mintime,
					skip_until_dynamic = self.is_dynamic(file) or skip_until_dynamic)
				
				if file in produced:
					if produced[file].rule > rule:
						continue
					elif produced[file].rule < rule:
						# prefer this rule, hence go on below
						pass
					else:
						if ignore_ambiguity:
							# ignore this job but don't throw error
							logger.warning("Rules {rule1} and {} are ambigous for file {}, using {rule1}.".format(rule, file, rule1=produced[file].rule))
							continue
						raise AmbiguousRuleException(file, produced[file], job, 
						                             lineno = self.lineno, 
						                             snakefile = self.snakefile)
				produced[file] = job
				
			except (ProtectedOutputException, MissingInputException, CyclicGraphException) as ex:
				exceptions[file].append(ex)
		
		missing_input = self._get_missing_files(set(input) - produced.keys())
		if missing_input:
			_ex = list()
			for file, exs in exceptions.items():
				if file in missing_input:
					_ex.extend(exs)
			raise MissingInputException(
				rule = self,
				files = set(missing_input) - exceptions.keys(), 
				include = _ex,
				lineno = self.lineno, 
				snakefile = self.snakefile
			)

		# collect the jobs that will actually run (including pseudo-jobs)
		todo = {job for job in produced.values() if job.needrun or job.pseudo}
		
		need_run, reason = self._need_run(forcethis or forceall, todo, input, output, parentmintime, requested_output, pseudo)
		
		if need_run:
			# enforce running jobs that created temporary files
			for file, job in produced.items():
				if not job.needrun and file.is_temp() and not file.exists():
					job.needrun = True
					job.reason = "Output file needed by rule {}: {}".format(self, file)
					todo.add(job)
		
		protected_output = self._get_protected_output(output) if need_run else None
		if protected_output:
			raise ProtectedOutputException(self, protected_output, 
			                               lineno = self.lineno, 
			                               snakefile = self.snakefile)
	
		for f in input:
			f.need()
		
		wildcards = Namedlist(fromdict = wildcards)

		job = Job(
			self.workflow,
			rule = self, 
			message = self.get_message(input, output, wildcards),
			reason = reason if give_reason else None,
			input = input,
			output = output,
			wildcards = wildcards,
			threads = self.threads,
			depends = todo,
			dryrun = dryrun,
			touch = touch,
			needrun = need_run,
			pseudo = pseudo,
			dynamic_output = [o for o in self.output if o in self.dynamic]
		)

		jobs[(output, self)] = job
		return job

	@staticmethod
	def _get_protected_output(output):
		return [o for o in output if os.path.exists(o) and not os.access(o, os.W_OK)]

	def check(self):
		"""
		Check if rule is well defined.
		"""
		if self.output and not self.has_run():
			raise RuleException("Rule {} defines output but does not have a \"run\" definition.".format(self.name), lineno = self.lineno, snakefile = self.snakefile)

	def _need_run(self, force, todo, input, output, parentmintime, requested_output, pseudo):
		""" Return True if rule needs to be run. """
		if self.has_run():
			if force:
				return True, "Forced rule execution."
			if not output:
				return True, ""

			# translate dynamic output to concrete output
			concrete_output = []
			if not pseudo:
				for i, f in enumerate(output):
					_f = self.output[i]
					if self.is_dynamic(_f):
						if f == requested_output:
							# remove requested output as it is only a placeholder
							requested_output = None
						concrete_output.extend([IOFile.create(d, temp=f.is_temp(), protected=f.is_protected()) for d, _ in listfiles(_f)])
					else:
						concrete_output.append(f)
			output = concrete_output
			
			output_mintime = IOFile.mintime(output) or parentmintime

			if self._has_missing_files(output, requested_output):
				return True, "Missing output files: {}".format(", ".join(self._get_missing_files(output)))
			if output_mintime == None: # TODO remove this since it cannot happen?
				return True, "Missing output files: {}".format(", ".join(self._get_missing_files(output)))

			if todo:
				todo_output = set()
				for job in todo:
					if not job.pseudo or job.needrun:
						todo_output.update(job.output)
				if todo_output:
					return True, "Updated input files: {}".format(", ".join(set(input) & set(todo_output)))


			newer = [i for i in input if os.path.exists(i) and i.is_newer(output_mintime)]
			if newer:
				return True, "Input files newer than output files: {}".format(", ".join(newer))
			
			return False, ""
		return False, ""

	def get_run(self):
		""" Return the run method. """
		#return self.workflow.get_run(self)
		return self.run_func

	def has_run(self):
		""" Return True if rule has a run method. """
		return self.run_func != None

	def get_message(self, input, output, wildcards, showmessage = True):
		"""
		Get the message that shall be printed upon rule execution.
		
		Arguments
		input       -- the input of the rule
		output      -- the output of the rule
		wildcards   -- the wildcards of the rule
		showmessage -- whether a user defined message shall be printed instead if existing  
		"""
		if self.message and showmessage:
			variables = dict(globals())
			variables.update(locals())
			msg = self.message.format(**variables)
			return msg

		def showtype(orig_iofiles, iofiles):
			for i, iofile in enumerate(iofiles):
				if orig_iofiles[i] in self.dynamic:
					f = str(orig_iofiles[i])
					f += " (dynamic)"
				else:
					f = str(iofile)
				if iofile.is_temp():
					f += " (temporary)"
				if iofile.is_protected():
					f += " (protected)"
				yield f

		msg = "rule " + self.name
		if input or output:
			msg += ":"
		if input:
			msg += "\n\tinput: {}".format(", ".join(showtype(self.input, input)))
		if output:
			msg += "\n\toutput: {}".format(", ".join(showtype(self.output, output)))
		return msg

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
		return bestmatch, bestmatch_output

	def clone(self):
		return Rule(self)
	
	@staticmethod
	def get_wildcard_len(wildcards):
		"""
		Return the length of the given wildcard values.
		
		Arguments
		wildcards -- a dict of wildcards
		"""
		return sum(map(len, wildcards.values()))

	def __lt__(self, rule):
		comp = self.workflow.get_ruleorder().compare(self.name, rule.name)
		return comp < 0

	def __gt__(self, rule):
		comp = self.workflow.get_ruleorder().compare(self.name, rule.name)
		return comp > 0

	def __repr__(self):
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
