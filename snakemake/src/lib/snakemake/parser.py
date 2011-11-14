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
		return type, name, self

class Rule(State):
	"""
	This state parses a snakemake rule.
	"""
	def __init__(self):
		self.name = None
		self.init = True
	def parse(self, type, name):
		if type == tokenize.NAME and self.init:
			self.init = False
			self.name = name
			return tokenize.STRING, 'Controller.get_instance().add_rule(Rule("{}"))'.format(name), self
		if type == tokenize.NAME and name == 'input':
			return tokenize.NEWLINE, '\n', Input(self)
		if type == tokenize.NAME and name == 'output':
			return tokenize.NEWLINE, '\n', Output(self)
		if type == tokenize.NAME and name == 'run':
			return tokenize.STRING, '\ndef {}(input, output, wildcards)'.format(self.name), Run(self)
		return None, None, self

class Input(State):
	"""
	This state parses the input for a snakemake rule
	"""
	def __init__(self, rule):
		self.rule = rule

	def parse(self, type, name):
		if type == tokenize.OP and name == ':':
			return tokenize.STRING, 'Controller.get_instance().last_rule().add_input([', self
		if type == tokenize.NEWLINE:
			return tokenize.STRING, '])', self.rule
		return type, name, self


class Output(State):
	"""
	This state parses the output for a snakemake rule
	"""
	def __init__(self, rule):
		self.rule = rule

	def parse(self, type, name):
		if type == tokenize.OP and name == ':':
			return tokenize.STRING, 'Controller.get_instance().last_rule().add_output([', self
		if type == tokenize.NEWLINE:
			return tokenize.STRING, '])', self.rule
		return type, name, self

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
	return tokenize.untokenize(translate(open(filepath).readline))
