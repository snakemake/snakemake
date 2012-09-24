# -*- coding: utf-8 -*-

import os, sys, re, stat, shutil, random
from itertools import product
from snakemake.exceptions import MissingOutputException, IOException
from snakemake.logging import logger

__author__ = "Johannes Köster"

class IOFile(str):
	"""
	A file that is either input or output of a rule.
	"""
	_register = dict()

	@classmethod
	def clear(cls):
		cls._register.clear()
	
	@classmethod
	def create(cls, file, temp = False, protected = False, origin = None):	
		if not isinstance(file, str) and not type(file).__name__ == "function":
			raise ValueError("Input and output files have to be specified as strings or functions that return a string given the used wildcards as single argument.")
		if origin is None:
			origin = file
		fid = (origin, file)
		if fid in cls._register:
			obj = cls._register[fid]
		else:
			obj = IOFile(file)
			cls._register[fid] = obj

		obj._temp = temp or obj._temp
		obj._protected = protected or obj._protected

		return obj

	@staticmethod
	def mintime(iofiles):
		existing = [f.mtime() for f in iofiles if os.path.exists(f)]
		if existing:
			return min(existing)
		return None
	
	def __init__(self, file):
		super().__init__(file)
		self._is_function = type(file).__name__ == "function"
		self._file = file
		self._regex = None
		self._needed = 0
		self._temp = False
		self._protected = False
				
	def get_file(self):
		if not self._is_function:
			return self._file
		else:
			raise ValueError("This IOFile is specified as a function and may not be used directly.")

	def is_temp(self):
		return self._temp

	def is_protected(self):
		return self._protected

	def need(self):
		self._needed += 1
	
	def used(self):
		self._needed -= 1
		if self._temp and self._needed == 0 and os.path.exists(self.get_file()):
			logger.warning("Deleting temporary file {}".format(self))
			os.remove(self.get_file())

	def exists(self):
		return os.path.exists(self.get_file())

	def mtime(self):
		return os.stat(self.get_file()).st_mtime
	
	def is_newer(self, time):
		return self.mtime() >= time
	
	def prepare(self):
		dir = os.path.dirname(self.get_file())
		if len(dir) > 0 and not os.path.exists(dir):
			os.makedirs(dir)
	
	def created(self, rulename, lineno, snakefile):
		if not os.path.exists(self.get_file()):
			raise MissingOutputException("Output file {} not produced by rule {}.".format(self.get_file(), rulename), lineno = lineno, snakefile = snakefile)
		if self._protected:
			logger.warning("Write protecting output file {}".format(self))
			mode = os.stat(self.get_file()).st_mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH
			if os.path.isdir(self.get_file()):
				for root, dirs, files in os.walk(self.get_file()):
					for d in dirs:
						os.chmod(os.path.join(self.get_file(), d), mode)
					for f in files:
						os.chmod(os.path.join(self.get_file(), f), mode)
			else:
				os.chmod(self.get_file(), mode)
	
	def remove(self):
		remove(self.get_file())
	
	def touch(self, rulename, lineno, snakefile):
		try:
			touch(self.get_file())
		except OSError as e:
			if e.errno == 2:
				raise MissingOutputException("Output file {} of rule {} shall be touched but does not exist.".format(self.get_file(), rulename), lineno = lineno, snakefile = snakefile)
			else:
				raise e


	def apply_wildcards(self, wildcards):
		f = self._file
		if self._is_function:
			f = self._file(Namedlist(fromdict = wildcards))
		return self.create(re.sub(_wildcard_regex, lambda match: '{}'.format(wildcards[match.group('name')]), f), protected = self._protected, temp = self._temp, origin = self)

	def fill_wildcards(self):
		f = self._file
		if self._is_function:
			raise ValueError("Cannot fill wildcards of function.")
		return self.create(re.sub(_wildcard_regex, lambda match: '{}'.format(random.randint(1, 1000000)), f), protected = self._protected, temp = self._temp)
		
	def get_wildcard_names(self):
		return set(match.group('name') for match in re.finditer(_wildcard_regex, self.get_file()))

	def regex(self):
		if not self._regex:
			# create a regular expression
			self._regex = re.compile(regex(self._file))
		return self._regex

	def match(self, target):
		match = re.match(self.regex(), target)
		if match and len(match.group()) == len(target):
			return match
		return None

_wildcard_regex = "\{\s*(?P<name>\w+?)(\s*,\s*(?P<constraint>[^\}]*))?\s*\}"

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

def touch(file):
	os.utime(self.get_file(), None)
	
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
	if not isinstance(filepatterns, list) or isinstance(filepatterns, tuple):
		filepatterns = [filepatterns]
	def flatten(wildcards):
		for wildcard, values in wildcards.items():
			if not isinstance(values, list) or isinstance(values, tuple):
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
	
	def __hash__(self):
		return hash(tuple(self))

	def __str__(self):
		return " ".join(self)

