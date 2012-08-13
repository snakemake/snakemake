# -*- coding: utf-8 -*-

import os, re
from snakemake.io import regex, Namedlist

__author__ = "Johannes Köster"

def linecount(filename):
	"""
	Return the number of lines of given file
	
	Arguments
	filename -- the path to the file
	"""
	with open(filename) as f:
		return sum(1 for l in f)

def listdir(dirname, pattern = None):
	"""
	Yield a tuple of filenames and full filepaths in the given directory.
	If pattern is specified, wildcard values are yielded as the third tuple item.

	Arguments
	dirname -- the path to the directory to list
	pattern -- an optional filepattern relative to dirname. Wildcards are specified in snakemake syntax, e.g. "{id}.txt"
	"""
	if pattern:
		pattern = re.compile(regex(pattern))
		for f in os.listdir(dirname):
			match = re.match(pattern, f)
			if match and len(match.group()) == len(f):
				wildcards = Namedlist(fromdict = match.groupdict())
				filepath = os.path.join(dirname, f)
				yield f, filepath, wildcards
	else:
		for f in os.listdir(dirname):
			filepath = os.path.join(dirname, f)
			yield f, filepath

def makedirs(dirnames):
	"""
	Recursively create the given directory or directories without reporting errors if they are present.
	"""
	if isinstance(dirnames, str):
		dirnames = [dirnames]
	for dirname in dirnames:
		if not os.path.exists(dirname):
			os.makedirs(dirname)
