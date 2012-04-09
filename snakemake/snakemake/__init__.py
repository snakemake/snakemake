# -*- coding: utf-8 -*-

import os, traceback, sys, csv
from snakemake.workflow import workflow
from snakemake.exceptions import print_exception
from snakemake.logging import logger, ColorizingStreamHandler

__author__ = "Johannes Köster"

def snakemake(snakefile, list = False, jobs = 1, directory = None, targets = None, dryrun = False, touch = False, forcethis = False, forceall = False, stats = None, give_reason = False, nocolor = False, cluster = None, standalone = False):
	"""
	Run snakemake on a given snakefile.
		
	Arguments
	snakefile         -- the snakefile.
	list              -- list rules.
	jobs              -- maximum number of parallel jobs (default: 1).
	directory         -- working directory (default: current directory).
	rule              -- execute this rule (default: first rule in snakefile).
	dryrun            -- print the rules that would be executed, but do not execute them.
	forcethis         -- force the selected rule to be executed
	forceall          -- force all rules to be executed
	time_measurements -- measure the running times of all rules
	"""

	if nocolor:
		ColorizingStreamHandler.nocolor = True

	def print_rules(log):
		log("Defined rules:")
		for rule in workflow.get_rules(): log(rule.name)
	
	olddir = os.getcwd()

	if directory:
		# change to the specified directory. This overrides eventually specified workdir in Snakefile
		workflow.set_workdir(directory)

	workflow.clear()

	if standalone:
		try:
			# set the process group
			os.setpgrp()
		except:
			# ignore: if it does not work we can still work without it
			pass

	try:
		workflow.include(snakefile, overwrite_first_rule = True)

		workflow.check_rules()

		if list:
			print_rules(logger.info)
			return 0
	
		workflow.set_cores(jobs)
		
		ret = 0
		if not targets: 
			ret = workflow.run_first_rule(dryrun = dryrun, touch = touch, forcethis = forcethis, forceall = forceall, give_reason = give_reason, cluster = cluster)
		else:
			ret = workflow.run_rules(targets, dryrun = dryrun, touch = touch, forcethis = forcethis, forceall = forceall, give_reason = give_reason, cluster = cluster)
		if ret == 0 and stats:
			stats = csv.writer(open(stats, "w"), delimiter = "\t")
			stats.writerow("rule minimum maximum sum mean".split())
			s = 0
			for measurement in workflow.get_runtimes():
				stats.writerow(measurement)
				s += measurement[3]
			stats.writerow([])
			stats.writerow(("Overall runtime", s))
		os.chdir(olddir)

		if standalone and ret == 1:
			try:
				# make sure ill behaving child processes are really killed (this will fail if snakemake is called programatically since it will kill the whole process)
				os.killpg(0, signal.SIGKILL)
			except:
				# ignore: if it does not work we can still work without it, but it may happen that some processes continue to run
				pass
		return ret
	except (Exception, BaseException) as ex:
		print_exception(ex, workflow.rowmaps)
		os.chdir(olddir)
		return 1
