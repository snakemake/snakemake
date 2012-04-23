# -*- coding: utf-8 -*-

import os, re
from snakemake.exceptions import MissingOutputException, IOException
from snakemake.logging import logger

__author__ = "Johannes Köster"

class IOFile(str):
	"""
	A file that is either input or output of a rule.
	"""
	wildcard_regex = "\{(?P<name>\w+?)(,(?P<constraint>.*))?\}"
	_register = dict()
	
	@classmethod
	def create(cls, file, temp = False, protected = False):	
		if not isinstance(file, str) and not type(file).__name__ == "function":
			raise ValueError("Input and output files have to be specified as strings or functions that return a string given the used wildcards as single argument.")
		if file in cls._register:
			obj = cls._register[file]
		else:
			obj = IOFile(file)
			cls._register[file] = obj

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
		self._is_function = type(file).__name__ == "function"
		self._file = file
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
		if self._temp and self._needed == 0:
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
			mode = os.stat(o).st_mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH
			if os.path.isdir(o):
				for root, dirs, files in os.walk(o):
					for d in dirs:
						os.chmod(os.path.join(o, d), mode)
					for f in files:
						os.chmod(os.path.join(o, f), mode)
			else:
				os.chmod(o, mode)
	
	def remove(self):
		if os.path.exists(self.get_file()):
			if os.path.isdir(self.get_file()):
				for root, dirs, files in os.walk(self.get_file()):
					for f in files:
						os.remove(f)
				for root, dirs, files in os.walk(self.get_file()):
					for d in dirs:
						os.rmdir(d)
			else:
				os.remove(self.get_file())
	
	def touch(self, rulename, lineno, snakefile):
		if not os.path.exists(self.get_file()):
			raise MissingOutputException("Output file {} of rule {} shall be touched but does not exist.".format(self.get_file(), rulename), lineno = lineno, snakefile = snakefile)
		os.utime(self.get_file(), None)

	def apply_wildcards(self, wildcards):
		f = self._file
		if self._is_function:
			f = self._file(wildcards)
		return self.create(re.sub(self.wildcard_regex, lambda match: '{}'.format(wildcards[match.group('name')]), f), protected = self._protected, temp = self._temp)
		
	def get_wildcard_names(self):
		return set(match.group('name') for match in re.finditer(self.wildcard_regex, self.get_file()))
	
	def regex(self):
		f = ""
		last = 0
		for match in re.finditer(self.wildcard_regex, self.get_file()):
			f += re.escape(self.get_file()[last:match.start()])
			f += "(?P<{}>{})".format(match.group("name"), match.group("constraint") if match.group("constraint") else ".+")
			last = match.end()
		f += re.escape(self.get_file()[last:])
		return f
	
	def __str__(self):
		return self.get_file()

	
class temp(str):
	pass

class protected(str):
	pass

