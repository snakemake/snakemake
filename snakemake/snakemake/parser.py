# -*- coding: utf-8 -*-

from tokenize import TokenError
from tokenize import *
import tokenize, collections, inspect

__author__ = "Johannes Köster"

class Tokens:
	""" Recorder for emitted tokens. """
	def __init__(self):
		self._row, self._col = 1, 0
		self._tokens = []
		self.rowmap = dict()
		
	def add(self, token, string = None, orig_token = None):
		""" Add a new token. Maybe a full TokenInfo object (first arg)
		or a pair of token type (e.g. NEWLINE) and string. """
		if string:
			type = token
			token = tokenize.TokenInfo(
				type = type,
				string = string,
				start = (self._row, self._col),
				end = (self._row, self._col + len(string)),
				line = '')
		else:
			# token is an original token that may have a wrong row
			if token.start[0] != self._row:
				token = Tokens._adjrow(token, self._row)
		
		self._tokens.append(token)
				
		if orig_token:
			self.rowmap[self._row] = orig_token.start[0]
						
		if token.type in (NEWLINE, NL):
			self._row += 1
			self._col = 0
		else:
			lines = token.string.split("\n")
			self._row += len(lines) - 1
			self._col += len(lines[-1]) + 1
		
		return self
	
	@staticmethod
	def _adjrow(token, row):
		""" Force the row of a token to be of the given value. """
		add = row - token.start[0]
		return token._replace(start = (token.start[0] + add, token.start[1]), 
			end = (token.end[0] + add, token.end[1]))
	
	def __iter__(self):
		return self._tokens.__iter__()

	def __str__(self):
		return " ".join(map(str, self))

class States:
	""" A finite automaton that translates snakemake tokens into python tokens. """
	def __init__(self, filename, rule_count = 0):
		self.state = self.python
		self.filename = filename
		self.main_states = dict(
			include = self.include,
			workdir = self.workdir,
			ruleorder = self.ruleorder,
			rule = self.rule,
			input = self.input,
			output = self.output,
			message = self.message,
			threads = self.threads,
			run = self.run,
			shell = self.shell)
		self.rule_params = set(["input", "output", "message", "threads", "run", "shell"])
		self.current_rule = None
		self.empty_rule = True
		self.tokens = Tokens()
		self._rule_count = rule_count
	
	def get_rule_count(self):
		return self._rule_count
	
	def __iter__(self):
		return self.tokens.__iter__()

	def _syntax_error(self, msg, token):
		""" Provide a convenient SyntaxError """
		return SyntaxError(msg, (self.filename, token.start[0], None, None))

	def python(self, token):
		""" The automaton state that handles ordinary python code. """
		if token.type == NAME and token.string in ('include', 'workdir', 'ruleorder', 'rule'):
			self.tokens.add(NEWLINE, '\n', orig_token = token)
			self.state = self.main_states[token.string]
		else:
			self.tokens.add(token, orig_token = token)
	
	def include(self, token):
		""" State that handles include definitions. """
		self._check_colon('include', token)
		self.state = self.include_path
		if self.empty_rule and self.current_rule:
			self.close_empty_rule(token)

	def include_path(self, token):
		""" State that translates the include path into a function call. """
		if token.type == STRING:
			self._func('include', (token.string,), token, obj = 'workflow')
			self.state = self.python
		else:
			raise self._syntax_error('Expected string after include keyword', token)

	def workdir(self, token):
		""" State that handles workdir definition. """
		self._check_colon('workdir', token)
		self.state = self.workdir_path
		if self.empty_rule and self.current_rule:
			self.close_empty_rule(token)

	def workdir_path(self, token):
		""" State that translates the workdir path into a function call. """
		if token.type == STRING:
			self._func('workdir', (token.string,), token, obj = 'workflow')
			self.state = self.python
		else:
			raise self._syntax_error('Expected string after workdir keyword', token)

	def ruleorder(self, token):
		""" State that handles ruleorder definitions. """
		self._check_colon('ruleorder', token)
		if self.empty_rule and self.current_rule:
			self.close_empty_rule(token)
		self._func_open('ruleorder', token, obj = 'workflow')
		self.state = self.ruleorder_order
	
	def ruleorder_order(self, token):
		if token.type == OP and token.string == '>':
			self.tokens.add(COMMA, ",", orig_token = token)
		elif token.type == NAME:
			self.tokens.add(STRING, self._stringify(token.string), orig_token = token)
		elif token.type == NEWLINE or token.type == NL or token.type == ENDMARKER:
			self._func_close(token)
			self.tokens.add(NEWLINE, '\n', orig_token = token)
			self.state = self.python
		elif not token.type in (INDENT, DEDENT, COMMENT):
			self._syntax_error('Expected a descending order of rule names, e.g. rule1 > rule2 > rule3 ...', token)

	def rule(self, token):
		""" State that handles rule definition. """
		if self.empty_rule and self.current_rule:
			self.close_empty_rule(token)

		self._rule_count += 1
		if self._is_colon(token):
			name = str(self._rule_count)
			self.state = self.rule_body
		elif token.type == NAME:
			name = token.string
			self.state = self.rule_colon
		else:
			raise self._syntax_error('Expected name or colon after rule keyword.', token)
		self.current_rule = name
		self.empty_rule = True
		self.tokens.add(AT, "@", token)
		self._func("rule", (self._stringify(name), str(token.start[0]), self._stringify(self.filename)), token, obj = 'workflow')
	
	def rule_colon(self, token):
		self._check_colon('rule', token)
		self.state = self.rule_body

	def rule_body(self, token):
		""" State that handles the rule body. """
		if token.type == NEWLINE or token.type == STRING:
			pass
		elif token.type == NAME and token.string in ("input", "output", "run", "shell", "threads", "message"):
			self.state = self.main_states[token.string]
		elif not token.type in (INDENT, DEDENT, COMMENT, NL):
			raise self._syntax_error('Expected one of the keywords "input", "output", "run", "shell", "threads" or "message" below rule definition', token)

	def input(self, token):
		""" State that handles input definition. """
		self.inoutput(token, 'input')

	def output(self, token):
		""" State that handles output definition. """
		self.inoutput(token, 'output')

	def inoutput(self, token, type):
		""" State that handles in- and output definition (depending on type). """
		self._check_colon(type, token)
		#self._func_open('set_{}'.format(type), token, obj = 'workflow')
		self.tokens.add(NEWLINE, "\n", token)\
		           .add(AT, "@", token)
		self._func_open(type, token, obj = 'workflow')
		           
		self.state = self.inoutput_paths
		
	def inoutput_paths(self, token):
		""" State that collects the arguments for in- or output definition """
		if token.type == NAME and token.string in self.main_states:
			self.state = self.main_states[token.string]
			self._func_close(token)
		elif token.type == ENDMARKER:
			self._func_close(token)
		else:
			self.tokens.add(token, orig_token = token)

	def message(self, token):
		""" State that handles message definition. """
		self._check_colon('message', token)
		self.tokens.add(NEWLINE, "\n", token)\
		           .add(AT, "@", token)
		self._func_open('message', token, obj = 'workflow')
		#self._func_open('set_message', token, obj = 'workflow')
		self.state = self.message_text

	def message_text(self, token):
		if token.type == STRING:
			self.tokens.add(token, orig_token = token)
			self.state = self.close_param
		elif not token.type in (INDENT, DEDENT, NEWLINE, NL):
			raise self._syntax_error('Expected string after message keyword.', token)
	
	def threads(self, token):
		""" State that handles definition of threads. """
		self._check_colon('thread', token)
		self.tokens.add(NEWLINE, "\n", token)\
		           .add(AT, "@", token)
		self._func_open('threads', token, obj = 'workflow')
		#self._func_open('set_threads', token, obj = 'workflow')
		self.state = self.threads_value
	
	def threads_value(self, token):
		if token.type == NUMBER:
			self.tokens.add(token, orig_token = token)
			self.state = self.close_param
		elif not token.type in (INDENT, DEDENT, NEWLINE, NL):
			raise self._syntax_error('Expected number after threads keyword.', token)
			

	def _run_def(self, token):
		self.tokens.add(NEWLINE, "\n", token)\
		           .add(AT, "@", token)\
		           .add(NAME, 'workflow', token)\
		           .add(DOT, '.', token)\
		           .add(NAME, 'run', token)\
		           .add(NEWLINE, '\n', token)
		self._func_def("__" + self.current_rule, ['input', 'output', 'wildcards', 'threads'], token)

	def run(self, token):
		""" State that creates a run function for the current rule. """
		self.empty_rule = False
		self._check_colon('run', token)
		self._run_def(token)
		self.state = self.run_newline

	def run_newline(self, token):
		if token.type == NEWLINE:
			self.tokens.add(token, orig_token = token)
		else:
			self.tokens.add(NEWLINE, '\n', orig_token = token)\
			           .add(INDENT, '\t', orig_token = token)\
			           .add(token, orig_token = token)
		self.state = self.python

	def run_body(self, token):
		""" State that collects the body of a rule's run function. """
		if token.type == NAME and token.string == 'rule':
			self.tokens.add(NEWLINE, '\n', orig_token = token)
			self.state = self.rule
		else:
			self.tokens.add(token, orig_token = token)

	def shell(self, token):
		""" State that creates a run function for the current rule, interpreting shell commands directly. """
		self.empty_rule = False
		self._check_colon('shell', token)
		self._run_def(token)
		
		self.tokens.add(NEWLINE, '\n', orig_token = token)\
		           .add(INDENT, '\t', orig_token = token)
		self._func_open('shell', token)
		self.state = self.shell_body

	def shell_body(self, token):
		""" State that collects the body of a rule's shell function. """
		if token.type == STRING:
			self.tokens.add(token, orig_token = token)
			self._func_close(token)
			self.state = self.python
		elif token.type in (COMMENT, NEWLINE, NL, INDENT, DEDENT):
			self.tokens.add(token, orig_token = token)
		else:
			raise self._syntax_error('Expected shell command in a string after shell keyword.', token)
	
	def close_param(self, token):
		""" State that closes a function invocation based definition of a rule parameter. """
		if token.type == NAME and token.string in self.main_states:
			self.state = self.main_states[token.string]
			self._func_close(token)
		elif token.type == ENDMARKER:
			self._func_close(token)

	def close_empty_rule(self, token):
		# close previous rule if empty
		self._run_def(token)
		self.tokens.add(NEWLINE, '\n', token)\
		           .add(INDENT, '\t', token)\
		           .add(NAME, 'pass', token)\
		           .add(NEWLINE, '\n', token)

	def _is_colon(self, token):
		return token.type == tokenize.OP and token.string == ':'

	def _check_colon(self, keyword, token):
		""" Check wether the token is a colon, else raise a syntax error 
		for the given keyword. """
		if not self._is_colon(token):
			raise self._syntax_error('Expected ":" after {} keyword'.format(keyword), token)

	def _func_def(self, name, args, orig_token):
		""" Generate tokens for a function definition with given name 
		and args. """
		self.tokens.add(NAME, 'def', orig_token = orig_token) \
				   .add(NAME, name, orig_token = orig_token) \
				   .add(LPAR, '(', orig_token = orig_token)
		for i in range(len(args)):
			self.tokens.add(NAME, args[i])
			if i < len(args) - 1: 
				self.tokens.add(COMMA, ',', orig_token = orig_token)
		self.tokens.add(RPAR, ')', orig_token = orig_token) \
				   .add(COLON, ':', orig_token = orig_token)

	def _func(self, name, args, orig_token, obj = None):
		""" Generate tokens for a function invocation with given name 
		and args. """
		if not isinstance(args, tuple):
			args = tuple(args)
		self._func_open(name, orig_token, obj = obj)
		for arg in args:
			self.tokens.add(STRING, arg, orig_token = orig_token) \
			           .add(COMMA, ',', orig_token = orig_token)
		self._func_close(orig_token)
		
	def _func_open(self, name, orig_token, obj = None):
		""" Generate tokens for opening a function invocation with 
		given name. """
		if obj != None:
			if obj:
				self.tokens.add(NAME, obj, orig_token = orig_token)
			self.tokens.add(DOT, '.', orig_token = orig_token)
		self.tokens.add(NAME, name, orig_token = orig_token) \
				   .add(LPAR, '(', orig_token = orig_token)

	def _func_close(self, orig_token):
		""" Generate tokens for closing a function invocation with 
		given name. """
		self.tokens.add(RPAR, ')', orig_token = orig_token) #\
		#		   .add(NEWLINE, '\n', orig_token = orig_token)

	@staticmethod
	def _stringify(tokenstring):
		""" Encapsulate a string into additional quotes. """
		return '"{}"'.format(tokenstring)			

def snakemake_to_python(tokens, filepath, rowmap = None, rule_count = 0):
	""" Translate snakemake tokens into python tokens using 
	a finite automaton. """
	states = States(filepath, rule_count = rule_count)
	try:
		for snakemake_token in tokens:
			states.state(snakemake_token)
	except TokenError as ex:
		raise SyntaxError(str(ex))
	python_tokens = (python_token for python_token in states if not python_token.type in (INDENT, DEDENT))
	if rowmap != None:
		rowmap.update(states.tokens.rowmap)
	return python_tokens, states.get_rule_count()

def compile_to_python(filepath, rule_count = 0):
	""" Compile a given Snakefile into python code. """
	with open(filepath) as snakefile:
		rowmap = dict()
		python_tokens, rule_count = snakemake_to_python(
				tokenize.generate_tokens(snakefile.readline), 
				filepath,
				rowmap = rowmap,
				rule_count = rule_count
		)
		compilation = tokenize.untokenize(python_tokens)
		return compilation, rowmap, rule_count
