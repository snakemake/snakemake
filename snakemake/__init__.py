# -*- coding: utf-8 -*-

import os
import argparse
from argparse import ArgumentError
import logging
import multiprocessing
import re
import sys

from snakemake.workflow import Workflow
from snakemake.exceptions import print_exception
from snakemake.logging import logger, init_logger

__author__ = "Johannes Köster"
__version__ = "2.4.5"


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
                    notemp=notemp,
                    cleanup_metadata=cleanup_metadata
                    )

    except (Exception, BaseException) as ex:
        print_exception(ex, workflow.linemaps)
    if workdir:
        os.chdir(olddir)
    if workflow.persistence:
        workflow.persistence.unlock()
    return success


def parse_resources(args):
    resources = dict()
    if args.resources is not None:
        valid = re.compile("[a-zA-Z_]\w*$")
        for res in args.resources:
            try:
                res, val = res.split("=")
            except ValueError:
                raise ValueError("Resources have to be defined as name=value pairs.")
            if not valid.match(res):
                raise ValueError("Resource definition must start with a valid identifier.")
            try:
                val = int(val)
            except ValueError:
                raise ValueError("Resource definiton must contain an integer after the identifier.")
            if res == "_cores":
                raise ValueError("Resource _cores is already defined internally. Use a different name.")
            resources[res] = val
    return resources


def main():
    parser = argparse.ArgumentParser(
        description="Snakemake is a Python based language and execution "
            "environment for GNU Make-like workflows.")

    parser.add_argument(
        "target", nargs="*", default=None,
        help="Targets to build. May be rules or files.")
    parser.add_argument(
        "--snakefile", "-s", metavar="FILE",
        default="Snakefile", help="The workflow definition in a snakefile.")
    parser.add_argument(
        "--cores", "--jobs", "-j", action="store", default=1,
        const=multiprocessing.cpu_count(), nargs="?", metavar="N", type=int,
        help=(
            "Use at most N cores in parallel (default: 1). "
            "If N is omitted, the limit is set to the number of "
            "available cores."))
    parser.add_argument(
        "--resources", "--res", nargs="*", metavar="NAME=INT",
        help=(
            "Define additional resources that shall constrain the scheduling "
            "analogously to threads (see above). A resource is defined as "
            "a name and an integer value. E.g. --resources gpu=1. Rules can "
            "use resources by defining the resource keyword, e.g. "
            "resources: gpu=1. If now two rules require 1 of the resource "
            "'gpu' they won't be run in parallel by the scheduler."))
    parser.add_argument(
        "--list", "-l", action="store_true",
        help="Show availiable rules in given snakefile.")
    parser.add_argument(
        "--directory", "-d", metavar="DIR", action="store",
        help=(
            "Specify working directory (relative paths in "
            "the snakefile will use this as their origin)."))
    parser.add_argument(
        "--dryrun", "-n", action="store_true",
        help="Do not execute anything.")
    parser.add_argument(
        "--printshellcmds", "-p", action="store_true",
        help="Print out the shell commands that will be executed.")
    parser.add_argument(
        "--dag", action="store_true",
        help="Do not execute anything and print the directed "
            "acyclic graph of jobs in the dot language. Recommended "
            "use on Unix systems: snakemake --dag | dot | display")
    parser.add_argument(
        "--ruledag", action="store_true",
        help="Do not execute anything and print the directed "
            "acyclic graph of rules in the dot language. This will be less "
            "crowded than above DAG of jobs, but also show less information. "
            "Use this if above option leads to a DAG that is too large. "
            "Recommended use on Unix systems: snakemake --ruledag | dot | display")
    parser.add_argument(
        "--summary", "-S", action="store_true",
        help="Print a summary of all files created by the workflow. The "
        "has the following columns: filename, modification time, "
        "rule version, status, plan.\n"
        "Thereby rule version contains the version"
        "the file was created with (see the version keyword of rules), and "
        "status denotes whether the file is missing, its input files are "
        "newer or if version or implementation of the rule changed since "
        "file creation. Finally the last column denotes whether the file "
        "will be updated or created during the next workflow execution.")
    parser.add_argument(
        "--touch", "-t", action="store_true",
        help=(
            "Touch output files (mark them up to date without really "
            "changing them) instead of running their commands. This is "
            "used to pretend that the rules were executed, in order to "
            "fool future invocations of snakemake. Fails if a file does "
            "not yet exist."))
    parser.add_argument(
        "--keep-going", "-k", action="store_true",
        help="Go on with independent jobs if a job fails.")
    parser.add_argument(
        "--force", "-f", action="store_true",
        help=(
            "Force the execution of the selected target or the first rule "
            "regardless of already created output."))
    parser.add_argument(
        "--forceall", "-F", action="store_true",
        help=(
            "Force the execution of the selected (or the first) rule and "
            "all rules it is dependent on regardless of already created "
            "output."))
    parser.add_argument(
        "--forcerun", "-R", nargs="+", metavar="TARGET",
        help=(
            "Force the re-execution or creation of the given rules or files."
            " Use this option if you changed a rule and want to have all its "
            "output in your workflow updated."))
    parser.add_argument(
        "--prioritize", "-P", nargs="+", metavar="TARGET",
        help=
            ("Tell the scheduler to assign creation of given targets "
            "(and all their dependencies) highest priority. (EXPERIMENTAL)"))
    parser.add_argument(
        "--allow-ambiguity", "-a", action="store_true",
        help=(
            "Don't check for ambiguous rules and simply use the first if "
            "several can produce the same file. This allows the user to "
            "prioritize rules by their order in the snakefile."))
    # TODO extend below description to explain the wildcards that can be used
    parser.add_argument(
        "--cluster", "-c", metavar="CMD",
        help=(
            "Execute snakemake rules with the given submit command, "
            "e.g. qsub. Snakemake compiles jobs into scripts that are "
            "submitted to the cluster with the given command, once all input "
            "files for a particular job are present (also see the argument below).\n"
            "The submit command can be decorated to make it aware of certain job properties (input, output, params, wildcards, log, threads and dependencies (see the argument below)), e.g.:\n"
            "$ snakemake --cluster 'qsub -pe threaded {threads}'."))
    parser.add_argument(
        "--immediate-submit", "--is", action="store_true",
        help=(
            "Immediately submit all jobs to the cluster instead of waiting "
            "for present input files. This will fail, unless you make "
            "the cluster aware of job dependencies, e.g. via:\n"
            "$ snakemake --cluster 'sbatch --dependency {dependencies}.\n"
            "Assuming that your submit script (here sbatch) outputs the generated job id to the first stdout line, {dependencies} will be filled with space separated job ids this job depends on."))
    parser.add_argument(
        "--jobscript", "--js", metavar="SCRIPT",
        help="Provide a custom job script for submission to the cluster. "
            "The default script resides as 'jobscript.sh' in the "
            "installation directory.")
    parser.add_argument(
        "--reason", "-r", action = "store_true",
        help="Print the reason for each executed rule.")
    parser.add_argument(
        "--stats", metavar="FILE",
        help="Write stats about Snakefile execution to the given file.")
    parser.add_argument(
        "--nocolor", action = "store_true",
        help="Do not use a colored output.")
    parser.add_argument(
        "--quiet", "-q", action = "store_true",
        help="Do not output any progress or rule information.")
    parser.add_argument(
        "--nolock", action="store_true",
        help="Do not lock the working directory")
    parser.add_argument(
        "--unlock", action="store_true",
        help="Remove a lock on the working directory.")
    parser.add_argument(
        "--cleanup-metadata", "--cm", nargs="*", metavar="FILE",
        help="Cleanup the metadata "
        "of given files. That means that snakemake removes any tracked "
        "version info, and any marks that files are incomplete.")
    parser.add_argument(
        "--rerun-incomplete", "--ri", action="store_true", help="Re-run all "
        "jobs the output of which is recognized as incomplete.")
    parser.add_argument(
        "--ignore-incomplete", "--ii", action="store_true", help="Ignore "
        "any incomplete jobs.")
    parser.add_argument(
        "--list-version-changes", "--lv", action="store_true",
        help="List all files that have been created with "
        "a different version (as determined by the version keyword).")
    parser.add_argument(
        "--list-code-changes", "--lc", action="store_true",
        help="List all files for which the rule body (run or shell) have changed "
        "in the Snakefile.")
    parser.add_argument(
        "--list-input-changes", "--li", action="store_true",
        help="List all files for which the defined input files have changed "
        "in the Snakefile.")
    parser.add_argument(
        "--list-params-changes", "--lp", action="store_true",
        help="List all files for which the defined params have changed "
        "in the Snakefile.")
    parser.add_argument(
        "--output-wait", "-w", type=int, default=3, metavar="SECONDS",
        help="Wait given seconds if an output file of a job is not present after "
        "the job finished. This helps if your filesystem "
        "suffers from latency.")
    parser.add_argument(
        "--notemp", "--nt", action="store_true",
        help="Ignore temp() declarations. This is useful when running only "
        "a part of the workflow, since temp() would lead to deletion of "
        "probably needed files by other parts of the workflow."
        )
    parser.add_argument(
        "--print-compilation", action="store_true",
        help="Print the python representation of the workflow.")
    parser.add_argument(
        "--debug", action="store_true", help="Print debugging output.")
    parser.add_argument(
        "--version", "-v", action="version", version=__version__)

    args = parser.parse_args()

    snakemakepath = os.path.realpath(__file__)

    try:
        resources = parse_resources(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        print("", file=sys.stderr)
        parser.print_help()
        return False

    success = snakemake(
            args.snakefile,
            listrules=args.list,
            cores=args.cores,
            resources=resources,
            workdir=args.directory,
            targets=args.target,
            dryrun=args.dryrun,
            printshellcmds=args.printshellcmds,
            printreason=args.reason,
            printdag=args.dag,
            printruledag=args.ruledag,
            touch=args.touch,
            forcetargets=args.force,
            forceall=args.forceall,
            forcerun=args.forcerun,
            prioritytargets=args.prioritize,
            stats=args.stats,
            nocolor=args.nocolor,
            quiet=args.quiet,
            keepgoing=args.keep_going,
            cluster=args.cluster,
            immediate_submit=args.immediate_submit,
            standalone=True,
            ignore_ambiguity=args.allow_ambiguity,
            snakemakepath=snakemakepath,
            lock=not args.nolock,
            unlock=args.unlock,
            cleanup_metadata=args.cleanup_metadata,
            force_incomplete=args.rerun_incomplete,
            ignore_incomplete=args.ignore_incomplete,
            list_version_changes=args.list_version_changes,
            list_code_changes=args.list_code_changes,
            list_input_changes=args.list_input_changes,
            list_params_changes=args.list_params_changes,
            summary=args.summary,
            print_compilation=args.print_compilation,
            debug=args.debug,
            jobscript=args.jobscript,
            notemp=args.notemp
            )
    return 0 if success else 1
