# -*- coding: utf-8 -*-

'''
Created on 13.11.2011

@author: Johannes Köster
'''

import tokenize

class State:
	"""
	An parser automaton state.
	"""
	def parse(self, type, name):
		"""
		Parse a token given by type and name and yield translated tokens.
		"""
		yield type, name, self

class Python(State):
	"""
	This state parses normal python code.
	"""
	def parse(self, type, name):
		if type == tokenize.NAME and name == 'rule':
			return tokenize.NEWLINE, '\n', Rule()
		if type == tokenize.NAME and name == 'workdir':
			return tokenize.NEWLINE, '\n', Workdir()
		return type, name, self

class Workdir(State):
	"""
	This state parses the definition of the working directory.
	"""
	def __init__(self):
		self.state = 0

	def parse(self, type, name):
		if type == tokenize.OP and name == ':':
			self.state = 1
			return tokenize.STRING, 'Controller.get_instance().set_workdir(', self
		elif self.state == 0:
			raise SyntaxError("Expected ':' after {} keyword".format(self.__class__.type))
		if type == tokenize.STRING and self.state == 1:
			self.state = 2
			return type, name, self
		if type == tokenize.NEWLINE:
			return tokenize.STRING, ')\n', Python()
		raise SyntaxError("Expected only one string after workdir keyword")

def is_keyword(name):
	return name in {"input", "output", "run", "rule"}

def handle_keywords(rule, name, prefix = ''):
	prefixtype = tokenize.NEWLINE
	if name == 'input':
		return prefixtype, prefix, Input(rule)
	if name == 'output':
		return prefixtype, prefix, Output(rule)
	if name == 'run':
		return tokenize.STRING, prefix + 'def {}(input, output, wildcards)'.format(rule.name), Run(rule)
	if name == 'rule':
		return prefixtype, prefix, Rule()
	raise SyntaxError("Unexpected keyword {} in rule".format(name))

class Rule(State):
	"""
	This state parses a snakemake rule.
	"""
	def __init__(self):
		self.name = None
		self.init = True
	def parse(self, type, name):
		if type == tokenize.NAME:
			if self.init:
				self.init = False
				self.name = name
				return tokenize.STRING, 'Controller.get_instance().add_rule(Rule("{}"))\n'.format(name), self
			return handle_keywords(self, name)
		return None, None, self

class Put(State):
	"""
	This state parses the input or output for a snakemake rule
	"""
	type = None
	def __init__(self, rule):
		self.rule = rule
		self.init = True

	def parse(self, type, name):
		if type == tokenize.OP and name == ':':
			self.init = False
			return tokenize.STRING, 'Controller.get_instance().last_rule().add_{}(['.format(self.__class__.type), self
		elif self.init:
			raise SyntaxError("Expected ':' after {} keyword".format(self.__class__.type))
		if is_keyword(name):
			return handle_keywords(self.rule, name, prefix = '])\n')
		if type == tokenize.ENDMARKER:
			return tokenize.STRING, '])\n', self
		return type, name, self

class Input(Put):
	type = 'input'

class Output(Put):
	type = 'output'

class Run(State):
	"""
	This state parses the run method for a snakemake rule
	"""
	def __init__(self, rule):
		self.rule = rule
		self.indent = 0

	def parse(self, type, name):
		if type == tokenize.INDENT:
			self.indent += 1
		if type == tokenize.DEDENT:
			self.indent -= 1
		if self.indent == -1:
			return tokenize.STRING, '', Python()
		return type, name, self

def translate(readline):
	"""
	Translate a snakemake DSL file to python code tokens
	"""
	state = Python()
	yield tokenize.NEWLINE, '\n' # fixes missing first token when untokenizing
	for type, name,_,_,_ in tokenize.generate_tokens(readline):
		type, name, state = state.parse(type, name)
		if type != None:
			yield type, name

def compile_to_python(filepath):
	"""
	Compile snakemake DSL to python code
	"""
	translated = list(translate(open(filepath).readline))
	compilation = tokenize.untokenize(translated)
	return compilation
