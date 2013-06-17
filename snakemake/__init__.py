# -*- coding: utf-8 -*-

import os

from snakemake.workflow import Workflow
from snakemake.exceptions import print_exception
from snakemake.logging import logger, init_logger

__author__ = "Johannes Köster"
__version__ = "2.4.3"


def snakemake(snakefile,
    listrules=False,
    cores=1,
    resources=None,
    workdir=None,
    targets=None,
    dryrun=False,
    touch=False,
    forcetargets=False,
    forceall=False,
    forcerun=None,
    prioritytargets=None,
    stats=None,
    printreason=False,
    printshellcmds=False,
    printdag=False,
    printruledag=False,
    nocolor=False,
    quiet=False,
    keepgoing=False,
    cluster=None,
    immediate_submit=False,
    standalone=False,
    ignore_ambiguity=False,
    snakemakepath=None,
    lock=True,
    unlock=False,
    cleanup_metadata=None,
    force_incomplete=False,
    ignore_incomplete=False,
    list_version_changes=False,
    list_code_changes=False,
    list_input_changes=False,
    list_params_changes=False,
    summary=False,
    output_wait=3,
    print_compilation=False,
    debug=False,
    notemp=False,
    jobscript=None):
    """
    Run snakemake on a given snakefile.
    Note: at the moment, this function is not thread-safe!

    Arguments
    snakefile         -- the snakefile.
    list              -- list rules.
    jobs              -- maximum number of parallel jobs (default: 1).
    directory         -- working directory (default: current directory).
    rule              -- execute this rule (default: first rule in snakefile).
    dryrun            -- print the rules that would be executed,
        but do not execute them.
    forcethis         -- force the selected rule to be executed
    forceall          -- force all rules to be executed
    time_measurements -- measure the running times of all rules
    lock              -- lock the working directory
    """

    init_logger(nocolor=nocolor, stdout=dryrun, debug=debug)

    if not os.path.exists(snakefile):
        logger.error("Error: Snakefile \"{}\" not present.".format(snakefile))
        return False

    if workdir:
        olddir = os.getcwd()
    workflow = Workflow(
        snakefile=snakefile, snakemakepath=snakemakepath,
        jobscript=jobscript)

    if standalone:
        try:
            # set the process group
            os.setpgrp()
        except:
            # ignore: if it does not work we can still work without it
            pass

    success = False
    try:
        workflow.include(snakefile, workdir=workdir,
            overwrite_first_rule=True, print_compilation=print_compilation)
        workflow.check()

        if not print_compilation:
            if listrules:
                workflow.list_rules()
            elif cleanup_metadata:
                workflow.cleanup_metadata(cleanup_metadata)
            else:
                success = workflow.execute(
                    targets=targets, dryrun=dryrun, touch=touch,
                    cores=cores, forcetargets=forcetargets,
                    forceall=forceall, forcerun=forcerun,
                    prioritytargets=prioritytargets, quiet=quiet,
                    keepgoing=keepgoing, printshellcmds=printshellcmds,
                    printreason=printreason, printruledag=printruledag,
                    printdag=printdag, cluster=cluster,
                    immediate_submit=immediate_submit,
                    ignore_ambiguity=ignore_ambiguity,
                    workdir=workdir, stats=stats,
                    force_incomplete=force_incomplete,
                    ignore_incomplete=ignore_incomplete,
                    list_version_changes=list_version_changes,
                    list_code_changes=list_code_changes,
                    list_input_changes=list_input_changes,
                    list_params_changes=list_params_changes,
                    summary=summary,
                    output_wait=output_wait,
                    nolock=not lock,
                    unlock=unlock,
                    resources=resources,
                    notemp=notemp
                    )

    except (Exception, BaseException) as ex:
        print_exception(ex, workflow.linemaps)
    if workdir:
        os.chdir(olddir)
    if workflow.persistence:
        workflow.persistence.unlock()
    return success
