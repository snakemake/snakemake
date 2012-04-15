# -*- coding: utf-8 -*-

import os, re
from snakemake.exceptions import MissingOutputException, IOException
from snakemake.logging import logger

__author__ = "Johannes Köster"

class IOFile(str):
	wildcard_regex = "\{(?P<name>\w+?)(,(?P<constraint>.*))?\}"
	_register = dict()
	
	@classmethod
	def create(cls, file, temp = False, protected = False):
		if not isinstance(file, str):
			print(file)
			raise ValueError("Input and output files have to be specified as strings.")

		obj = None
		if file in cls._register:
			obj = cls._register[file]
		else:
			obj = IOFile(file)
			cls._register[file] = obj

		if temp:	
			obj._temp = temp
		if protected:
			obj._protected = protected
		return obj

	@staticmethod
	def mintime(iofiles):
		existing = [f.mtime() for f in iofiles if os.path.exists(f)]
		if existing:
			return min(existing)
		return None
	
	def __init__(self, file):
		self._file = file
		self._needed = 0
		self._temp = False
		self._protected = False

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
			os.remove(self._file)

	def exists(self):
		return os.path.exists(self._file)

	def mtime(self):
		return os.stat(self._file).st_mtime
	
	def is_newer(self, time):
		return self.mtime() >= time
	
	def prepare(self):
		dir = os.path.dirname(self._file)
		if len(dir) > 0 and not os.path.exists(dir):
			os.makedirs(dir)
	
	def created(self, rulename, lineno, snakefile):
		if not os.path.exists(self._file):
			raise MissingOutputException("Output file {} not produced by rule {}.".format(self._file, rulename), lineno = lineno, snakefile = snakefile)
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
		if os.path.exists(self._file):
			if os.path.isdir(self._file):
				for root, dirs, files in os.walk(self._file):
					for f in files:
						os.remove(f)
				for root, dirs, files in os.walk(self._file):
					for d in dirs:
						os.rmdir(d)
			else:
				os.remove(self._file)
	
	def touch(self, rulename, lineno, snakefile):
		if not os.path.exists(self._file):
			raise MissingOutputException("Output file {} of rule {} shall be touched but does not exist.".format(self._file, rulename), lineno = lineno, snakefile = snakefile)
		os.utime(self._file, None)

	def apply_wildcards(self, wildcards):
		return self.create(re.sub(self.wildcard_regex, lambda match: '{}'.format(wildcards[match.group('name')]), self._file), protected = self._protected, temp = self._temp)
		
	def get_wildcard_names(self):
		return set(match.group('name') for match in re.finditer(self.wildcard_regex, self._file))
	
	def regex(self):
		f = re.sub("\.", "\.", self._file)
		return re.sub(self.wildcard_regex, 
			lambda match: '(?P<{}>{})'.format(match.group('name'), match.group('constraint') if match.group('constraint') else ".+"), f)
	
	def __str__(self):
		return self._file

class temp(str):
	pass

class protected(str):
	pass
