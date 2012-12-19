# -*- coding: utf-8 -*-

import os, io, re, fnmatch, mimetypes, base64, inspect, textwrap, tempfile, subprocess, shutil, mimetypes, datetime
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

def report(text, path, stylesheet = os.path.join(os.path.dirname(__file__), "report.css"), defaultenc="utf8", template=None, **files):
	outmime, _ = mimetypes.guess_type(path)
	if outmime != "text/html":
		raise ValueError("Path to report output has to be an HTML file.")
	from docutils.core import publish_file
	metadata = textwrap.dedent("""

	.. raw:: html

	   <div id="metadata">{}</div>
	
	""").format(datetime.date.today().isoformat())
	text = format(textwrap.dedent(text), stepout=2) + metadata
	attachments = []
	for name, file in files.items():
		mime, encoding = mimetypes.guess_type(file)
		if mime is None:
			mime = ""
		if encoding is None:
			encoding = defaultenc
		with open(file, "rb") as f:
			data = base64.b64encode(f.read())
		attachments.append('.. _{}: data:{};charset={};base64,{}'.format(name, mime, encoding, data.decode()))
	text += "\n\n" + "\n\n".join(attachments)
	overrides = dict()
	if template is not None:
		overrides["template"] = template
	if stylesheet is not None:
		overrides["stylesheet_path"] = stylesheet
	html = open(path, "w")
	publish_file(source=io.StringIO(text), destination=html, writer_name="html", settings_overrides=overrides)

def R(code):
	import rpy2.robjects as robjects
	robjects.r(format(textwrap.dedent(code), stepout=2))

def format(string, *args, stepout = 1, **kwargs):
	class SequenceFormatter:
		def __init__(self, sequence):
			self._sequence = sequence

		def __getitem__(self, i):
			return self._sequence[i]

		def __str__(self):
			return " ".join(self._sequence)
		
	frame = inspect.currentframe().f_back
	while stepout > 1:
		if not frame.f_back:
			break
		frame = frame.f_back
		stepout -= 1
	
	variables = dict(frame.f_globals)
	# add local variables from calling rule/function
	variables.update(frame.f_locals)
	variables.update(kwargs)
	strmethods = list()
	for key, value in list(variables.items()):
		if type(value) in (list, tuple, set, frozenset):
			variables[key] = SequenceFormatter(value)
	try:
		return string.format(*args, **variables)
	except KeyError as ex:
		raise NameError("The name {} is unknown in this context.".format(str(ex)))
