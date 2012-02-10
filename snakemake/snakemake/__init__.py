import os, logging, traceback, sys, csv
from snakemake.workflow import workflow
from snakemake.parser import compile_to_python
from snakemake.exceptions import print_exception

def snakemake(snakefile, list = False, jobs = 1, directory = None, targets = None, dryrun = False, forcethis = False, forceall = False, stats = None):
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
	
	def print_rules(log):
		log("Defined rules:")
		for rule in workflow.get_rules(): log(rule.name)
	
	code, rowmap = compile_to_python(snakefile)

	if directory:
		# change to the specified directory. This overrides eventually specified workdir in Snakefile
		workflow.set_workdir(directory)    

	workflow.clear()

	try:
		workflow.execdsl(code, rowmap)

		workflow.check_rules()

		if list:
			print_rules(logging.info)
			return 0
	
		workflow.setup_pool(jobs)
		
		if not targets: 
			workflow.run_first_rule(dryrun = dryrun, forcethis = forcethis, forceall = forceall)
		else:
			workflow.run_rules(targets, dryrun = dryrun, forcethis = forcethis, forceall = forceall)
		if stats:
			stats = csv.writer(open(stats, "w"), delimiter = "\t")
			stats.writerow("rule minimum maximum sum mean".split())
			s = 0
			for measurement in workflow.get_runtimes():
				stats.writerow(measurement)
				s += measurement[3]
			stats.writerow([])
			stats.writerow(("Overall runtime", s))
	except (Exception, BaseException) as ex:
		print_exception(ex, rowmap)
		return 1
	return 0
