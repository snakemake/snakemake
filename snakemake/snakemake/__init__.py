import os, sys, logging
from snakemake.workflow import workflow, RuleException
from snakemake.parser import compile_to_python

def snakemake(snakefile, list = False, jobs = 1, directory = None, target = None, dryrun = False, forcethis = False, forceall = False):
	"""
	Run snakemake on a given snakefile.
		
	Arguments
	snakefile -- the snakefile.
	list -- list rules.
	jobs -- maximum number of parallel jobs (default: 1).
	directory -- working directory (default: current directory).
	rule -- execute this rule (default: first rule in snakefile).
	dryrun -- print the rules that would be executed, but do not execute them.
	"""
	def print_rules(log):
		log("Defined rules:")
		for rule in workflow.get_rules(): log(rule.name)
	
	code = compile_to_python(snakefile)

	if directory:
		# change to the specified directory. This overrides eventually specified workdir in Snakefile
		workflow.set_workdir(directory)    

	workflow.execdsl(code)

	workflow.check_rules()

	if list:
		print_rules(logging.info)
		return 0
	
	workflow.setup_pool(jobs)
	try:
		if not target: 
			workflow.run_first_rule(dryrun = dryrun, forcethis = forcethis, forceall = forceall)
		else:
			if workflow.is_rule(target): 
				workflow.run_rule(target, dryrun = dryrun, forcethis = forcethis, forceall = forceall)
			else:
				workflow.produce_file(target, dryrun = dryrun, forcethis = forcethis, forceall = forceall)
	except RuleException as ex:
		logging.error(str(ex))
		return 1
	return 0
