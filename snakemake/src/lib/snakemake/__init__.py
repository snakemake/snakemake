import os, sys, logging
from snakemake.workflow import Controller, RuleException
from snakemake.parser import compile_to_python

def snakemake(snakefile, list=False, jobs=1, directory=None, rule=None, dryrun=False, force=False):
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
    Controller.jobs = jobs
    controller = Controller.get_instance()
    
    def print_rules(log):
        log("Defined rules:")
        for rule in controller.get_rules(): log(rule.name)
    
    code = compile_to_python(snakefile)
    
    if directory:
        if os.path.exists(directory): os.chdir(directory)
        else:
            logging.error("Error: Defined working directory does not exist.")
            return 1
    
    controller.execdsl(code)
    
    if list:
        print_rules(logging.info)
        return 0
    
    #controller.setup_dag()
    try:
        if not rule: controller.apply_first_rule(dryrun=dryrun, force=force)
        elif controller.is_rule(rule): controller.apply_rule(rule, dryrun=dryrun, force=force)
        else:
            logging.error("Error: Rule {} does not exist.\n".format(rule))
            print_rules(logging.error)
            return 1
    except RuleException as ex:
        logging.error(str(ex))
