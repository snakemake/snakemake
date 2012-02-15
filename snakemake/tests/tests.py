import sys
import os
from os.path import join
from subprocess import call
from tempfile import mkdtemp
import hashlib
from snakemake import snakemake

__author__ = "Tobias Marschall, Marcel Martin"


def dpath(path):
	"""get path to a data file (relative to the directory this
	test lives in)"""
	return join(os.path.dirname(__file__), path)


def md5sum(filename):
	data = open(filename, 'rb').read()
	return hashlib.md5(data).hexdigest()


def run(path, **params):
	"""
	Test the Snakefile in path.
	There must be a Snakefile in the path and a subdirectory named
	expected-results.
	"""
	results_dir = join(path, 'expected-results')
	snakefile = join(path, 'Snakefile')
	assert os.path.exists(snakefile)
	assert os.path.exists(results_dir) and os.path.isdir(results_dir), \
		'{} does not exist'.format(results_dir)
	tmpdir = mkdtemp()
	try:
		call('cp `find {} -maxdepth 1 -type f` {}'.format(path, tmpdir), shell=True)
		# TODO
		# The snakemake call changes the current working directory, so
		# we need to remember it and restore it below.
		olddir = os.getcwd()
		exitcode = snakemake(snakefile, directory=tmpdir, stats = os.path.join(path, "stats.txt"),**params)
		os.chdir(olddir)
		assert exitcode == 0, "exit code is not zero, but {}".format(exitcode)
		for resultfile in os.listdir(results_dir):
			if not os.path.isfile(resultfile):
				continue # skip .svn dirs etc.
			targetfile = join(tmpdir, resultfile)
			expectedfile = join(results_dir, resultfile)
			assert os.path.exists(targetfile), 'expected file "{}" not produced'.format(resultfile)
			assert md5sum(targetfile) == md5sum(expectedfile), 'wrong result produced for file "{}"'.format(resultfile)
	finally:
		call(['rm', '-rf', tmpdir])


def test01():
	run(dpath("test01"))


def test02():
	run(dpath("test02"))


def test03():
	run(dpath("test03"), targets=['test.out'])


def test04():
	run(dpath("test04"), targets=['test.out'])

def test06():
	run(dpath("test06"), targets=['test.bla.out'])

def test07():
	run(dpath("test07"), targets=['test.out', 'test2.out'])

def test08():
	run(dpath("test08"))
