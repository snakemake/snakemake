# -*- coding: utf-8 -*-

import sys, traceback
from collections import defaultdict
from snakemake.logging import logger
from tokenize import TokenError

__author__ = "Johannes Köster"


def format_error(ex, lineno, rowmaps = None, snakefile = None):
	msg = str(ex)
	if rowmaps and snakefile:
		lineno = rowmaps[snakefile][lineno]
		if isinstance(ex, SyntaxError):
			msg = ex.msg
	location = " in line {} of {}".format(lineno, snakefile) if lineno and snakefile else ""
	return '{}{}{}'.format(ex.__class__.__name__, location, ":\n" + msg if msg else ".")

def get_exception_origin(ex, rowmaps):
	for file, lineno, _, _ in reversed(traceback.extract_tb(ex.__traceback__)):
		#traceback.print_tb(ex.__traceback__)
		if file in rowmaps:
			return lineno, file

def print_exception(ex, rowmaps):
	"""
	Print an error message for a given exception.

	Arguments
	ex -- the exception
	rowmaps -- a dict of a dict that maps for each snakefile the compiled lines to source code lines in the snakefile.
	"""
	origin = get_exception_origin(ex, rowmaps)
	if not origin is None:
		lineno, file = origin
		logger.critical(format_error(ex, lineno, rowmaps = rowmaps, snakefile = file))
		return
	if isinstance(ex, SyntaxError):
		logger.critical(format_error(ex, ex.lineno, rowmaps = rowmaps, snakefile = ex.filename))
	elif isinstance(ex, RuleException):
		for e in ex._include + [ex]:
			if not e.omit:
				logger.critical(format_error(e, e.lineno, rowmaps = rowmaps, snakefile = e.filename))
	elif isinstance(ex, KeyboardInterrupt):
		logger.warning("Cancelling snakemake on user request.")
	else:
		traceback.print_tb(ex.__traceback__)
		logger.critical(ex)

class RuleException(Exception):
	"""
	Base class for exception occuring withing the execution or definition of rules.
	"""
	def __init__(self, message = None, include = list(), lineno = None, snakefile = None):
		"""
		Creates a new instance of RuleException.

		Arguments
		message -- the exception message
		include -- iterable of other exceptions to be included
		lineno -- the line the exception originates
		snakefile -- the file the exception originates
		"""
		super(RuleException, self).__init__(message)
		self._include = set()
		for ex in include:
			self._include.add(ex)
			self._include.update(ex._include)
		self._include = list(self._include)
		self.lineno = lineno
		self.filename = snakefile
		self.omit = not message

class MissingOutputException(RuleException):
	pass
		
class IOException(RuleException):
	def __init__(self, prefix, rule, files, include = list(), lineno = None, snakefile = None):
		message = "{} for rule {}:\n{}".format(prefix, rule, "\n".join(files)) if files else ""
		super().__init__(message = message, include = include, lineno = lineno, snakefile = snakefile)

class MissingInputException(IOException):
	def __init__(self, rule, files, include = list(), lineno = None, snakefile = None):
		super().__init__("Missing input files", rule, files, include, lineno = lineno, snakefile = snakefile)

class ProtectedOutputException(IOException):
	def __init__(self, rule, files, include = list(), lineno = None, snakefile = None):
		super().__init__("Write-protected output files", rule, files, include, lineno = lineno, snakefile = snakefile)

class AmbiguousRuleException(RuleException):
	def __init__(self, filename, job1, job2, lineno = None, snakefile = None):
		super().__init__("Rules {} and {} are ambiguous for the file {}. Advice: Use --dag to see the two alternative DAGs.".format(job1.rule, job2.rule, filename), lineno = lineno, snakefile = snakefile)
		self.job1, self.job2 = job1, job2
	
class CyclicGraphException(RuleException):
	def __init__(self, rule, lineno = None, snakefile = None):
		super().__init__("Cyclic dependency on rule {}.".format(rule), lineno = lineno, snakefile = snakefile)
		
class MissingRuleException(RuleException):
	def __init__(self, file, lineno = None, snakefile = None):
		super().__init__("No rule to produce {}.".format(file), lineno = lineno, snakefile = snakefile)

class UnknownRuleException(RuleException):
	def __init__(self, name, lineno = None, snakefile = None):
		super().__init__("There is no rule named {}.".format(name), lineno = lineno, snakefile = snakefile)
		
class NoRulesException(RuleException):
	def __init__(self, lineno = None, snakefile = None):
		super().__init__("There has to be at least one rule.", lineno = lineno, snakefile = snakefile)

class IOFileException(RuleException):
	def __init__(self, msg, lineno = None, snakefile = None):
		super().__init__(msg, lineno = lineno, snakefile = snakefile)

class ClusterJobException(RuleException):
	def __init__(self, job):
		super().__init__("Error executing rule {} on cluster. For detailed error see the cluster log.".format(job.rule.name), lineno = job.rule.lineno, snakefile = job.rule.snakefile)

class CreateRuleException(RuleException):
	pass

class TerminatedException(Exception):
	pass
