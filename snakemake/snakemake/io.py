# -*- coding: utf-8 -*-

import os, sys, re, stat, shutil, random
from itertools import product
from snakemake.exceptions import MissingOutputException, IOException
from snakemake.logging import logger

__author__ = "Johannes Köster"

def IOFile(file, rule=None):
	f = _IOFile(file)
	f.rule = rule
	return f

class _IOFile(str):
	"""
	A file that is either input or output of a rule.
	"""
	
	def __new__(cls, file):
		obj = str.__new__(cls, file)
		obj._is_function = type(file).__name__ == "function"
		obj._file = file
		obj.rule = None
		obj._regex = None
		return obj
				
	@property
	def file(self):
		if not self._is_function:
			return self._file
		else:
			raise ValueError("This IOFile is specified as a function and may not be used directly.")

	@property
	def exists(self):
		return os.path.exists(self.file)
	
	@property
	def protected(self):
		return self.exists and not os.access(self.file, os.W_OK)

	@property
	def mtime(self):
		return os.stat(self.file).st_mtime
	
	def is_newer(self, time):
		return self.mtime > time
	
	def prepare(self):
		path_until_wildcard = re.split(_wildcard_regex, self.file)[0]
		dir = os.path.dirname(path_until_wildcard)
		if len(dir) > 0 and not os.path.exists(dir):
			try:
				os.makedirs(dir)
			except OSError as e:
				# ignore Errno 17 "File exists" (may happen due to multiprocessing)
				if e.errno != 17:
					raise e
	
	def protect(self):
		logger.warning("Write protecting output file {}".format(self))
		mode = os.stat(self.file).st_mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH
		if os.path.isdir(self.file):
			for root, dirs, files in os.walk(self.file):
				for d in dirs:
					os.chmod(os.path.join(self.file, d), mode)
				for f in files:
					os.chmod(os.path.join(self.file, f), mode)
		else:
			os.chmod(self.file, mode)
	
	def remove(self):
		remove(self.file)
	
	def touch(self):
		try:
			os.utime(self.file, None)
		except OSError as e:
			if e.errno == 2:
				raise MissingOutputException("Output file {} of rule {} shall be touched but does not exist.".format(self.file, self.rule.name), lineno = self.rule.lineno, snakefile = self.rule.snakefile)
			else:
				raise e
	
	def apply_wildcards(self, wildcards):
		f = self._file
		if self._is_function:
			f = self._file(Namedlist(fromdict = wildcards))
		return IOFile(re.sub(_wildcard_regex, lambda match: '{}'.format(wildcards[match.group('name')]), f), rule=self.rule)
	
	def fill_wildcards(self):
		f = self._file
		if self._is_function:
			raise ValueError("Cannot fill wildcards of function.")
		return IOFile(re.sub(_wildcard_regex, lambda match: "0", f), rule=self.rule)
		
	def get_wildcard_names(self):
		return set(match.group('name') for match in re.finditer(_wildcard_regex, self.file))

	def regex(self):
		if not self._regex:
			# create a regular expression
			self._regex = re.compile(regex(self.file))
		return self._regex

	def match(self, target):
		match = re.match(self.regex(), target)
		if match and len(match.group()) == len(target):
			return match
		return None
	
	def __eq__(self, other):
		f = other._file if isinstance(other, _IOFile) else other
		return self._file == f
	
	def __hash__(self):
		return self._file.__hash__()

_wildcard_regex = re.compile("\{\s*(?P<name>\w+?)(\s*,\s*(?P<constraint>[^\}]*))?\s*\}")

def remove(file):
	if os.path.exists(file):
		if os.path.isdir(file):
			try:
				os.removedirs(file)
			except OSError:
				# ignore non empty directories
				pass
		else:
			os.remove(file)
	
def regex(filepattern):
	f = ""
	last = 0
	for match in re.finditer(_wildcard_regex, filepattern):
		f += re.escape(filepattern[last:match.start()])
		f += "(?P<{}>{})".format(match.group("name"), match.group("constraint") if match.group("constraint") else ".+")
		last = match.end()
	f += re.escape(filepattern[last:])
	return f

class temp(str):
	""" A flag for an input or output file that shall be removed after usage. """
	pass

class temporary(temp):
	""" An alias for temp. """
	pass

class protected(str):
	""" A flag for a file that shall be write protected after creation. """
	pass

class dynamic(str):
	""" A flag for a file that shall be dynamic, i.e. the multiplicity (and wildcard values) will be expanded after a certain rule has been run """
	def __new__(cls, file):
		matches = list(re.finditer(_wildcard_regex, file))
		if len(matches) != 1:
			raise SyntaxError("Dynamic files need exactly one wildcard.")
		if matches[0].group("constraint"):
			raise SyntaxError("The wildcard in dynamic files cannot be constrained.")
		obj = str.__new__(cls, file)
		return obj
	pass


def expand(*args, **wildcards):
	""" 
	Expand wildcards in given filepatterns. 

	Arguments
	*args -- first arg: filepatterns as list or one single filepattern, second arg (optional): a function to combine wildcard values (itertools.product per default)
	**wildcards -- the wildcards as keyword arguments with their values as lists
	"""
	filepatterns = args[0]
	if len(args) == 1:
		combinator = product
	elif len(args) == 2:
		combinator = args[1]
	if isinstance(filepatterns, str):
		filepatterns = [filepatterns]
	def flatten(wildcards):
		for wildcard, values in wildcards.items():
			if isinstance(values, str):
				values = [values]
			yield [(wildcard, value) for value in values]
	expanded = list()
	for comb in combinator(*flatten(wildcards)):
		comb = dict(comb)
		for filepattern in filepatterns:
			expanded.append(filepattern.format(**comb))
	return expanded
	
class Namedlist(list):
	"""
	A list that additionally provides functions to name items. Further,
	it is hashable, however the hash does not consider the item names.
	"""
	def __init__(self, toclone = None, fromdict = None):
		"""
		Create the object.
		
		Arguments
		toclone  -- another Namedlist that shall be cloned
		fromdict -- a dict that shall be converted to a Namedlist (keys become names) 
		"""
		super(Namedlist, self).__init__()
		self._names = dict()
		
		if toclone:
			self.extend(toclone)
			if isinstance(toclone, Namedlist):
				self.take_names(toclone.get_names())
		if fromdict:
			for key, item in fromdict.items():
				self.append(item)
				self.add_name(key)

	def add_name(self, name):
		"""
		Add a name to the last item.
		
		Arguments
		name -- a name
		"""
		self.set_name(name, len(self) - 1)
	
	def set_name(self, name, index):
		"""
		Set the name of an item.
		
		Arguments
		name  -- a name
		index -- the item index
		"""
		self._names[name] = index
		setattr(self, name, self[index])
			
	def get_names(self):
		"""
		Get the defined names as (name, index) pairs.
		"""
		for name, index in self._names.items():
			yield name, index
	
	def take_names(self, names):
		"""
		Take over the given names.
		
		Arguments
		names -- the given names as (name, index) pairs
		"""
		for name, index in names:
			self.set_name(name, index)

	def items(self):
		for name, index in self._names.items():
			yield name, self[index]

	def insert_items(self, index, items):
		self[index:index+1] = items
		for name, i in self._names.items():
			if i > index:
				self._names[name] += len(items) - 1
	
	def __hash__(self):
		return hash(tuple(self))

	def __str__(self):
		return " ".join(self)


