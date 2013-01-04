# -*- coding: utf-8 -*-

import os, traceback, sys, csv
from snakemake.workflow import Workflow
from snakemake.exceptions import print_exception
from snakemake.logging import logger, init_logger

__author__ = "Johannes Köster"
__version__ = "2.1.1"

def snakemake(snakefile, 
              listrules = False, 
              cores = 1, 
              workdir = None, 
              targets = None, 
              dryrun = False, 
              touch = False, 
              forcetargets = False, 
              forceall = False, 
              forcerules = None, 
              prioritytargets = None,
              stats = None, 
              printreason = False, 
              printshellcmds = False, 
              printdag = False, 
              nocolor = False, 
              quiet = False, 
              keepgoing = False,
              cluster = None, 
              standalone = False,
              ignore_ambiguity = False, 
              snakemakepath = None):
	"""
	Run snakemake on a given snakefile.
	Note: at the moment, this function is not thread-safe!
		
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
	init_logger(nocolor=nocolor, stdout=dryrun)

	if not os.path.exists(snakefile):
		logger.error("Error: Snakefile \"{}\" not present.".format(snakefile))
		return False
	
	if workdir:
		olddir = os.getcwd()
	workflow = Workflow(snakemakepath = snakemakepath)

	if standalone:
		try:
			# set the process group
			os.setpgrp()
		except:
			# ignore: if it does not work we can still work without it
			pass

	success = False
	try:
		workflow.include(snakefile, workdir=workdir, overwrite_first_rule=True)

		if listrules:
			workflow.list_rules()
			return True
		
		success = workflow.execute(targets=targets, dryrun=dryrun, touch=touch, 
		                           cores=cores, forcetargets=forcetargets, 
		                           forceall=forceall, forcerules=forcerules, 
		                           prioritytargets=prioritytargets, quiet=quiet, keepgoing=keepgoing,
		                           printshellcmds=printshellcmds, printreason=printreason, 
		                           printdag=printdag, cluster=cluster, 
		                           ignore_ambiguity=ignore_ambiguity, 
		                           workdir=workdir, stats=stats)
		
	except (Exception, BaseException) as ex:
		print_exception(ex, workflow.linemaps)
	if workdir:
		os.chdir(olddir)
	return success
