# -*- coding: utf-8 -*-

import os, re, fnmatch
from itertools import chain
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

def listfiles(pattern):
	"""
	Yield a tuple of existing filepaths for the given pattern.
	If pattern is specified, wildcard values are yielded as the third tuple item.

	Arguments
	pattern -- a filepattern. Wildcards are specified in snakemake syntax, e.g. "{id}.txt"
	"""
	first_wildcard = re.search("{[^{]", pattern)
	if first_wildcard:
		dirname = os.path.dirname(pattern[:first_wildcard.start()])
		if not dirname:
			dirname = "."
	else:
		dirname = os.path.dirname(pattern)
	pattern = re.compile(regex(pattern))
	for dirpath, dirnames, filenames in os.walk(dirname):
		for f in chain(filenames, dirnames):
			if dirpath != ".":
				f = os.path.join(dirpath, f)
			match = re.match(pattern, f)
			if match and len(match.group()) == len(f):
				wildcards = Namedlist(fromdict = match.groupdict())
				yield f, wildcards

def makedirs(dirnames):
	"""
	Recursively create the given directory or directories without reporting errors if they are present.
	"""
	if isinstance(dirnames, str):
		dirnames = [dirnames]
	for dirname in dirnames:
		if not os.path.exists(dirname):
			os.makedirs(dirname)
