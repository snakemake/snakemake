from token import *
import tokenize, collections, inspect

class States:
	''' A finite automaton that tranlates snakemake tokens into python tokens. '''
	def __init__(self, filename):
		self.state = self.python
		self.filename = filename
		self.main_states = dict(
			workdir = self.workdir,
			rule = self.rule,
			input = self.input,
			output = self.output,
			run = self.run)
		self.current_rule = None

	def _syntax_error(self, msg, token):
		''' Provide a convenient SyntaxError '''
		return SyntaxError(msg, (self.filename, token.start[0], None, None))

	def python(self, token):
		''' The automaton state that handles ordinary python code. '''
		if token.type == NAME:
			if token.string in ('workdir', 'rule'):
				self.state = self.main_states[token.string]
			else:
				yield token
		else:
			yield token

	def workdir(self, token):
		''' State that handles workdir definition. '''
		self._check_colon('workdir', token)
		self.state = self.workdir_path

	def workdir_path(self, token):
		''' State that translates the workdir path into a function call. '''
		if token.type == STRING:
			for t in States._func(token, '_set_workdir', token.string): yield t
			self.state = self.python
		else:
			raise self._syntax_error('Expected string after workdir keyword', token)

	def rule(self, token):
		''' State that handles rule definition. '''
		if token.type == NAME:
			self.current_rule = token.string
			for t in States._func(token, '_add_rule', States._stringify(token.string)): yield t
			self.state = self.rule_colon
		else:
			raise self._syntax_error('Expected name after rule keyword.', token)
	
	def rule_colon(self, token):
		self._check_colon('rule', token)
		self.state = self.rule_body

	def rule_body(self, token):
		''' State that handles the rule body. '''
		if token.type in (INDENT, DEDENT, NEWLINE):
			yield token
		elif token.type == NAME and token.string in self.main_states:
			self.state = self.main_states[token.string]
		else:
			print(token)
			raise self._syntax_error('Expected one of the keywords "input", "output", "run" or "rule" below rule definition', token)

	def input(self, token):
		''' State that handles input definition. '''
		self.inoutput(token, 'input')

	def output(self, token):
		''' State that handles output definition. '''
		self.inoutput(token, 'output')

	def inoutput(self, token, type):
		''' State that handles in- and output definition (depending on type). '''
		self._check_colon(type, token)
		for t in States._func_open(token, '_set_{}'.format(type)): yield t
		self.state = self.inoutput_paths
		
	def inoutput_paths(self, token):
		''' State that collects the arguments for in- or output definition '''
		if token.type == NAME and token.string in self.main_states:
			self.state = self.main_states[token.string]
			yield States._func_close(token)
		else:
			yield token

	def run(self, token):
		''' State that creates a run function for the current rule. '''
		self._check_colon(type, token)
		for t in States._func_def(token, self.current_rule,
			['input', 'output', 'wildcards']): yield t
		self.state = self.run_body

	def run_body(self, token):
		''' State that collects the body of a rule's run function. '''
		if token.type == NAME and token.string in self.main_states:
			self.state = self.main_states[token.string]
		else:
			yield token

	def _check_colon(self, keyword, token):
		''' Check wether the token is a colon, else raise a syntax error 
		for the given keyword. '''
		if not (token.type == tokenize.OP and token.string == ':'):
			raise self._syntax_error('Expected ":" after workdir keyword', token)

	@staticmethod
	def _func_def(token, name, args):
		''' Generate tokens for a function definition with given name 
		and args. '''
		yield token._replace(type = NAME, string = 'def')
		yield token._replace(type = NAME, string = name)
		yield token._replace(type = LPAR, string = '(')
		for i in range(len(args)):
			yield token._replace(type = NAME, string = arg)
			if i < len(args) - 1: yield token._replace(type = COMMA, string = ',')
		yield token._replace(type = RPAR, string = ')')
		yield token._replace(type = COLON, string = ':')

	@staticmethod
	def _func(token, name, arg):
		''' Generate tokens for a function invocation with given name 
		and args. '''
		token = States._fix_col(token)
		for t in States._func_open(token, name): yield t
		yield token._replace(type=STRING, string=arg)
		yield States._func_close(token)
		
	@staticmethod
	def _func_open(token, name):
		''' Generate tokens for opening a function invocation with 
		given name. '''
		yield token._replace(type = NAME, string = name)
		yield token._replace(type = LPAR, string = '(')

	@staticmethod
	def _func_close(token):
		''' Generate tokens for closing a function invocation with 
		given name. '''
		return token._replace(type = RPAR, string = ')')
	
	@staticmethod
	def _fix_col(token):
		return token._replace(start = (token.start[0], 0))
	
	@staticmethod
	def _stringify(tokenstring):
		return '"{}"'.format(tokenstring)


def snakemake_to_python(tokens, filepath):
	states = States(filepath)
	for snakemake_token in tokens:
		python_tokens = states.state(snakemake_token)
		if inspect.isgenerator(python_tokens):
			for python_token in python_tokens:
				yield python_token

def compile_to_python(filepath):
	with open(filepath) as snakefile:
		snakemake_tokens = list(snakemake_to_python(
			tokenize.generate_tokens(snakefile.readline), filepath))
		#print(snakemake_tokens)
		#print(len(snakemake_tokens))
		compilation = tokenize.untokenize(snakemake_tokens)
		print(compilation)
		return compilation
