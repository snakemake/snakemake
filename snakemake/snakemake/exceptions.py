import sys, traceback
from collections import defaultdict

def format_error(ex, lineno, rowmaps = None, snakefile = None):
	if rowmaps and snakefile:
		lineno = rowmaps[snakefile][lineno]
	return '{} in line {} of {}:\n{}'.format(ex.__class__.__name__, lineno, snakefile, str(ex))

def print_exception(ex, rowmaps):
	for file, lineno, _, _ in traceback.extract_tb(ex.__traceback__):
		if file in rowmaps:
			print(format_error(ex, lineno, rowmaps = rowmaps, snakefile = file), file = sys.stderr)
			return
	if isinstance(ex, SyntaxError):
		print(format_error(ex, ex.lineno, rowmaps = rowmaps, snakefile = ex.filename), file = sys.stderr)
	elif isinstance(ex, RuleException):
		for e in ex._include + [ex]:
			if not e.omit:
				print(format_error(e, e.lineno, rowmaps = rowmaps, snakefile = e.filename), file = sys.stderr)
	else:
		traceback.print_tb(ex.__traceback__)
		print(ex, file=sys.stderr)

class RuleException(Exception):
	def __init__(self, message = None, include = list(), lineno = None, snakefile = None):
		super(RuleException, self).__init__(message)
		self._include = set()
		for ex in include:
			self._include.add(ex)
			self._include.update(ex._include)
		self._include = list(self._include)
		self.lineno = lineno
		self.filename = snakefile
		self.omit = len(message) == 0

class MissingOutputException(RuleException):
	pass
		
class IOException(RuleException):
	def __init__(self, prefix, rule, files, include = list(), lineno = None, snakefile = None):
		message = "{} for rule {}:\n{}".format(prefix, rule, "\n".join(files)) if files else ""
		super(IOException, self).__init__(message = message, include = include, lineno = lineno, snakefile = snakefile)

class MissingInputException(IOException):
	def __init__(self, rule, files, include = list(), lineno = None, snakefile = None):
		super(MissingInputException, self).__init__("Missing input files", rule, files, include, lineno = lineno, snakefile = snakefile)

class ProtectedOutputException(IOException):
	def __init__(self, rule, files, include = list(), lineno = None, snakefile = None):
		super(ProtectedOutputException, self).__init__("Protected output files", rule, files, include, lineno = lineno, snakefile = snakefile)

class AmbiguousRuleException(RuleException):
	def __init__(self, rule1, rule2, lineno = None, snakefile = None):
		super(AmbiguousRuleException, self).__init__("Rules {} and {} are ambiguous.".format(rule1, rule2), lineno = lineno, snakefile = snakefile)

class CyclicGraphException(RuleException):
	def __init__(self, rule, lineno = None, snakefile = None):
		super(CyclicGraphException, self).__init__("Cyclic dependency on rule {}.".format(rule), lineno = lineno, snakefile = snakefile)
		
class MissingRuleException(RuleException):
	def __init__(self, file, lineno = None, snakefile = None):
		super(MissingRuleException, self).__init__("No rule to produce {}.".format(file), lineno = lineno, snakefile = snakefile)

class UnknownRuleException(RuleException):
	def __init__(self, name, lineno = None, snakefile = None):
		super(UnknownRuleException, self).__init__("There is no rule named {}.".format(name), lineno = lineno, snakefile = snakefile)
		
class NoRulesException(RuleException):
	def __init__(self, lineno = None, snakefile = None):
		super(NoRulesException, self).__init__("There has to be at least one rule.", lineno = lineno, snakefile = snakefile)

class CreateRuleException(RuleException):
	pass
