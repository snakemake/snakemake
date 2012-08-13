# -*- coding: utf-8 -*-

import os, re
from snakemake.io import regex, Namedlist

__author__ = "Johannes Köster"

def linecount(filename):
	with open(filename) as f:
		return sum(1 for l in f)

def listdir(dirname, pattern = None):
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
	if isinstance(dirnames, str):
		dirnames = [dirnames]
	for dirname in dirnames:
		if not os.path.exists(dirname):
			os.makedirs(dirname)
