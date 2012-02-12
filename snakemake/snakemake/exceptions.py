import sys, traceback
from collections import defaultdict

def print_exception(ex, rowmap):
	traceback.print_tb(ex.__traceback__)

	for file, lineno, _, _ in traceback.extract_tb(ex.__traceback__):
			if file == "<string>":
				print("Error in line {} of Snakefile:\n{}".format(rowmap[lineno], str(ex)), file = sys.stderr)
				return
	if not isinstance(ex, RuleException):
		traceback.print_tb(ex.__traceback__)
	print(ex, file=sys.stderr)

class RuleException(Exception):
	def __init__(self, message = None, include = list()):
		super(RuleException, self).__init__(message)
		self._include = set(ex for ex in include)
	
	def __str__(self):
		messages = set(str(ex) for ex in self._include)
		message = super(RuleException, self).__str__()
		if message: messages.add(message)
		return "\n".join(messages)

class MissingOutputException(RuleException):
	pass
		
class IOException(RuleException):
	def __init__(self, prefix, rule, files, include = list()):
		message = "{} for rule {}:\n{}".format(prefix, rule, "\n".join(files)) if files else ""
		super(IOException, self).__init__(message = message, include = include)

class MissingInputException(IOException):
	def __init__(self, rule, files, include = list()):
		super(MissingInputException, self).__init__("Missing input files", rule, files, include)

class ProtectedOutputException(IOException):
	def __init__(self, rule, files, include = list()):
		super(ProtectedOutputException, self).__init__("Protected output files", rule, files, include)

class AmbiguousRuleException(RuleException):
	def __init__(self, rule1, rule2):
		super(AmbiguousRuleException, self).__init__("Ambiguous rules: {} and {}.".format(rule1, rule2))

class CyclicGraphException(RuleException):
	def __init__(self, rule):
		super(CyclicGraphException, self).__init__("Cyclic dependency on rule {}.".format(rule))
		
class MissingRuleException(RuleException):
	def __init__(self, file):
		super(MissingRuleException, self).__init__("No rule to produce {}.".format(file))

class CreateRuleException(RuleException):
	pass
