#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from subprocess import call,Popen,PIPE
from tempfile import mkdtemp

__author__ = "Tobias Marschall"

usage = """{} <snakemake-executable>"""

def md5sum(filename):
	f = Popen(['md5sum', filename], stdout=PIPE)
	return f.stdout.readline().split()[0]

def main():
	if len(sys.argv) != 2:
		print(usage.format(sys.argv[0]), file=sys.stderr)
		sys.exit(1)
	executable = sys.argv[1]
	n = 0
	success = 0
	for f in os.listdir('.'):
		if not os.path.isdir(f): continue
		if not f.startswith('test'): continue
		if not os.path.exists(f+'/Snakefile'):
			print('Warning: {}/Snakefile does not exist, skipping directory'.format(f), file=sys.stderr)
			continue
		results_dir = '{}/expected-results'.format(f)
		if not os.path.exists(results_dir) or not os.path.isdir(results_dir):
			print('Warning: {} does not exist, skipping directory'.format(results_dir), file=sys.stderr)
			continue
		n += 1
		print('RUNNING {}'.format(f))
		tmpdir = mkdtemp()
		try:
			call('cp `find {} -maxdepth 1 -type f` {}'.format(f,tmpdir), shell=True)
			additional_params = ''
			if os.path.exists(f+'/params'):
				additional_params = open(f+'/params').readline().strip()
			exitcode = call('{0} -d {1} -s {1}/Snakefile {2} > /dev/null 2>&1'.format(executable,tmpdir,additional_params), shell=True)
			fail = False
			if exitcode != 0:
				print(' - FAILED with exit code', exitcode)
				continue
			fail = False
			for resultfile in os.listdir(results_dir):
				if not os.path.isfile(resultfile): continue
				targetfile = '{}/{}'.format(tmpdir,resultfile)
				expectedfile = '{}/{}'.format(results_dir,resultfile)
				if not os.path.exists(targetfile):
					print(' - ERROR: expected file "{}" not produced'.format(resultfile))
					fail = True
				if md5sum(targetfile) != md5sum(expectedfile):
					print(' - ERROR: wrong result produced for file "{}"'.format(resultfile))
					fail = True
			if not fail:
				print(' - SUCCESS')
				success += 1
		finally:
			call(['rm', '-rf', tmpdir])
	print('SUMMARY: {} of {} tests successful'.format(success, n))

if __name__ == '__main__':
	sys.exit(main())
