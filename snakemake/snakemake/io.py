import os, re
from snakemake.exceptions import MissingOutputException

class IOFile(str):
	wildcard_regex = "\{(?P<name>\w+?)(,(?P<constraint>.*))?\}"
	_register = dict()
	
	@classmethod
	def create(cls, file, temp = False):
		obj = None
		if file in cls._register:
			obj = cls._register[file]
		else:
			obj = IOFile(file)
			cls._register[file] = obj 
		if obj._temp == None:
			obj._temp = temp
		return obj
	
	def __init__(self, file, protected = False):
		self._file = file
		self._needed = 0
		self._temp = None
		self._protected = protected
		
	def need(self):
		self._needed += 1
	
	def used(self):
		self._needed -= 1
		if self._temp:
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
	
	def touch(self, lineno, snakefile):
		if not os.path.exists(self._file):
			raise MissingOutputException("Output file {} shall be touched but does not exist.".format(self._file), lineno = lineno, snakefile = snakefile)
		os.utime(self._file, None)

	def apply_wildcards(self, wildcards):
		return self.create(re.sub(self.wildcard_regex, lambda match: '{}'.format(wildcards[match.group('name')]), self._file))
		#return self.create(f.format(**wildcards))
		
	def get_wildcard_names(self):
		return set(match.group('name') for match in re.finditer(self.wildcard_regex, self._file))
	
	def regex(self):
		f = re.sub("\.", "\.", self._file)
		return re.sub(self.wildcard_regex, 
			lambda match: '(?P<{}>{})'.format(match.group('name'), match.group('constraint') if match.group('constraint') else ".+"), f)
	
	def __str__(self):
		return self._file

def temp(file):
	return IOFile.create(file, temp = True)

def protected(file):
	return IOFile.create(file, protected = True)
