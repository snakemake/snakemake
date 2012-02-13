import sys, traceback
from collections import defaultdict

def format_error(ex, lineno, rowmap = None):
	if rowmap:
		lineno = rowmap[lineno]
	return "Error in line {} of Snakefile:\n{}".format(lineno, str(ex))

def print_exception(ex, rowmap):
	for file, lineno, _, _ in traceback.extract_tb(ex.__traceback__):
		if file == "<string>":
			print(format_error(ex, lineno, rowmap = rowmap), file = sys.stderr)
			return
	if isinstance(ex, SyntaxError):
		print(format_error(ex, ex.lineno, rowmap = rowmap), file = sys.stderr)
		return
	if isinstance(ex, RuleException):
		if ex.lineno:
			print(format_error(ex, ex.lineno), file = sys.stderr)
			return
	
	if not isinstance(ex, RuleException):
		traceback.print_tb(ex.__traceback__)
	print(ex, file=sys.stderr)

class RuleException(Exception):
	def __init__(self, message = None, include = list(), lineno = None):
		super(RuleException, self).__init__(message)
		self._include = set(ex for ex in include)
		self.lineno = lineno
	
	def __str__(self):
		messages = set(str(ex) for ex in self._include)
		message = super(RuleException, self).__str__()
		if message: messages.add(message)
		return "\n".join(messages)

class MissingOutputException(RuleException):
	pass
		
class IOException(RuleException):
	def __init__(self, prefix, rule, files, include = list(), lineno = None):
		message = "{} for rule {} in line {}:\n{}".format(prefix, rule, lineno, "\n".join(files)) if files else ""
		super(IOException, self).__init__(message = message, include = include)

class MissingInputException(IOException):
	def __init__(self, rule, files, include = list(), lineno = None):
		super(MissingInputException, self).__init__("Missing input files", rule, files, include, lineno = lineno)

class ProtectedOutputException(IOException):
	def __init__(self, rule, files, include = list(), lineno = None):
		super(ProtectedOutputException, self).__init__("Protected output files", rule, files, include, lineno = lineno)

class AmbiguousRuleException(RuleException):
	def __init__(self, rule1, rule2, lineno = None):
		super(AmbiguousRuleException, self).__init__("Ambiguous rules: {} and {}.".format(rule1, rule2), lineno = lineno)

class CyclicGraphException(RuleException):
	def __init__(self, rule, lineno = None):
		super(CyclicGraphException, self).__init__("Cyclic dependency on rule {}.".format(rule), lineno = lineno)
		
class MissingRuleException(RuleException):
	def __init__(self, file, lineno = None):
		super(MissingRuleException, self).__init__("No rule to produce {}.".format(file), lineno = lineno)

class UnknownRuleException(RuleException):
	def __init__(self, name, lineno = None):
		super(UnknownRuleException, self).__init__("There is no rule named {}.".format(name), lineno = lineno)
		
class NoRulesException(RuleException):
	def __init__(self, lineno = None):
		super(NoRulesException, self).__init__("There has to be at least one rule.", lineno = lineno)

class CreateRuleException(RuleException):
	pass
