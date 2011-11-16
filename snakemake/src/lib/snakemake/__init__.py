import os, sys
from snakemake.workflow import Controller
from snakemake.parser import compile_to_python

def snakemake(snakefile, list=False, jobs=1, directory=None, rule=None, dryrun=False):
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
    
    def print_rules(file = sys.stdout):
        print("Defined rules:", file=file)
        for rule in controller.get_rules(): print(rule.name, file=file)
    
    code = compile_to_python(snakefile)
    
    if directory:
        if os.path.exists(directory): os.chdir(directory)
        else:
            print("Error: Defined working directory does not exist.", file = sys.stderr)
            return 1
    
    controller.execdsl(code)
    
    if list:
        print_rules()
        return 0
    
    controller.setup_dag()
    
    if not rule: controller.apply_first_rule(dryrun=dryrun)
    elif controller.is_rule(rule): controller.apply_rule(rule, dryrun=dryrun)
    else:
        print("Error: Rule {} does not exist.".format(rule), file=sys.stderr)
        print_rules(file=sys.stderr)
        return 1