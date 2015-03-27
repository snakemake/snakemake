# -*- coding: utf-8 -*-

import os
import subprocess
import glob
import argparse
from argparse import ArgumentError
import logging as _logging
import multiprocessing
import re
import sys
import inspect
import threading
import webbrowser
from functools import partial


from snakemake.workflow import Workflow
from snakemake.exceptions import print_exception
from snakemake.logging import setup_logger, logger
from snakemake.version import __version__
from snakemake.io import load_configfile
from snakemake.shell import shell


__author__ = "Johannes Köster"


def snakemake(snakefile,
    listrules=False,
    list_target_rules=False,
    cores=1,
    nodes=1,
    resources=dict(),
    config=dict(),
    configfile=None,
    config_args=None,
    workdir=None,
    targets=None,
    dryrun=False,
    touch=False,
    forcetargets=False,
    forceall=False,
    forcerun=[],
    prioritytargets=[],
    stats=None,
    printreason=False,
    printshellcmds=False,
    printdag=False,
    printrulegraph=False,
    printd3dag=False,
    nocolor=False,
    quiet=False,
    keepgoing=False,
    cluster=None,
    drmaa=None,
    jobname="snakejob.{rulename}.{jobid}.sh",
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
    list_resources=False,
    summary=False,
    detailed_summary=False,
    latency_wait=3,
    benchmark_repeats=1,
    wait_for_files=None,
    print_compilation=False,
    debug=False,
    notemp=False,
    nodeps=False,
    keep_target_files=False,
    allowed_rules=None,
    jobscript=None,
    timestamp=False,
    greedyness=None,
    overwrite_shellcmd=None,
    updated_files=None,
    log_handler=None,
    keep_logger=False):
    """Run snakemake on a given snakefile.

    This function provides access to the whole snakemake functionality. It is not thread-safe.

    Args:
        snakefile (str):            the path to the snakefile
        listrules (bool):           list rules (default False)
        list_target_rules (bool):   list target rules (default False)
        cores (int):                the number of provided cores (ignored when using cluster support) (default 1)
        nodes (int):                the number of provided cluster nodes (ignored without cluster support) (default 1)
        resources (dict):           provided resources, a dictionary assigning integers to resource names, e.g. {gpu=1, io=5} (default {})
        config (dict):              override values for workflow config
        workdir (str):              path to working directory (default None)
        targets (list):             list of targets, e.g. rule or file names (default None)
        dryrun (bool):              only dry-run the workflow (default False)
        touch (bool):               only touch all output files if present (default False)
        forcetargets (bool):        force given targets to be re-created (default False)
        forceall (bool):            force all output files to be re-created (default False)
        forcerun (list):             list of files and rules that shall be re-created/re-executed (default [])
        prioritytargets (list):     list of targets that shall be run with maximum priority (default [])
        stats (str):                path to file that shall contain stats about the workflow execution (default None)
        printreason (bool):         print the reason for the execution of each job (default false)
        printshellcmds (bool):      print the shell command of each job (default False)
        printdag (bool):            print the dag in the graphviz dot language (default False)
        printrulegraph (bool):      print the graph of rules in the graphviz dot language (default False)
        printd3dag (bool):          print a D3.js compatible JSON representation of the DAG (default False)
        nocolor (bool):             do not print colored output (default False)
        quiet (bool):               do not print any default job information (default False)
        keepgoing (bool):           keep goind upon errors (default False)
        cluster (str):              submission command of a cluster or batch system to use, e.g. qsub (default None)
        drmaa (str):                if not None use DRMAA for cluster support, str specifies native args passed to the cluster when submitting a job
        jobname (str):              naming scheme for cluster job scripts (default "snakejob.{rulename}.{jobid}.sh")
        immediate_submit (bool):    immediately submit all cluster jobs, regardless of dependencies (default False)
        standalone (bool):          kill all processes very rudely in case of failure (do not use this if you use this API) (default False)
        ignore_ambiguity (bool):    ignore ambiguous rules and always take the first possible one (default False)
        snakemakepath (str):        path to the snakemake executable (default None)
        lock (bool):                lock the working directory when executing the workflow (default True)
        unlock (bool):              just unlock the working directory (default False)
        cleanup_metadata (bool):    just cleanup metadata of output files (default False)
        force_incomplete (bool):    force the re-creation of incomplete files (default False)
        ignore_incomplete (bool):   ignore incomplete files (default False)
        list_version_changes (bool): list output files with changed rule version (default False)
        list_code_changes (bool):   list output files with changed rule code (default False)
        list_input_changes (bool):  list output files with changed input files (default False)
        list_params_changes (bool): list output files with changed params (default False)
        summary (bool):             list summary of all output files and their status (default False)
        latency_wait (int):         how many seconds to wait for an output file to appear after the execution of a job, e.g. to handle filesystem latency (default 3)
        benchmark_repeats (int):    number of repeated runs of a job if declared for benchmarking (default 1)
        wait_for_files (list):      wait for given files to be present before executing the workflow
        list_resources (bool):      list resources used in the workflow (default False)
        summary (bool):             list summary of all output files and their status (default False). If no option  is specified a basic summary will be ouput. If 'detailed' is added as an option e.g --summary detailed, extra info about the input and shell commands will be included
        detailed_summary (bool):    list summary of all input and output files and their status (default False)
        print_compilation (bool):   print the compilation of the snakefile (default False)
        debug (bool):               show additional debug output (default False)
        notemp (bool):              ignore temp file flags, e.g. do not delete output files marked as temp after use (default False)
        nodeps (bool):              ignore dependencies (default False)
        keep_target_files (bool):   Do not adjust the paths of given target files relative to the working directory.
        allowed_rules (set):        Restrict allowed rules to the given set. If None or empty, all rules are used.
        jobscript (str):            path to a custom shell script template for cluster jobs (default None)
        timestamp (bool):           print time stamps in front of any output (default False)
        greedyness (float):         set the greedyness of scheduling. This value between 0 and 1 determines how careful jobs are selected for execution. The default value (0.5 if prioritytargets are used, 1.0 else) provides the best speed and still acceptable scheduling quality.
        overwrite_shellcmd (str):   a shell command that shall be executed instead of those given in the workflow. This is for debugging purposes only.
        updated_files(list):        a list that will be filled with the files that are updated or created during the workflow execution
        log_handler (function):      redirect snakemake output to this custom log handler, a function that takes a log message dictionary (see below) as its only argument (default None). The log message dictionary for the log handler has to following entries:

            :level:
                the log level ("info", "error", "debug", "progress", "job_info")

            :level="info", "error" or "debug":
                :msg:
                    the log message
            :level="progress":
                :done:
                    number of already executed jobs

                :total:
                    number of total jobs

            :level="job_info":
                :input:
                    list of input files of a job

                :output:
                    list of output files of a job

                :log:
                    path to log file of a job

                :local:
                    whether a job is executed locally (i.e. ignoring cluster)

                :msg:
                    the job message

                :reason:
                    the job reason

                :priority:
                    the job priority

                :threads:
                    the threads of the job


    Returns:
        bool:   True if workflow execution was successful.

    """

    if updated_files is None:
        updated_files = list()

    if cluster:
        cores = sys.maxsize
    else:
        nodes = sys.maxsize

    if not keep_logger:
        setup_logger(handler=log_handler, quiet=quiet, printreason=printreason, printshellcmds=printshellcmds, nocolor=nocolor, stdout=dryrun, debug=debug, timestamp=timestamp)

    if greedyness is None:
         greedyness = 0.5 if prioritytargets else 1.0
    else:
        if not (0 <= greedyness <= 1.0):
            logger.error("Error: greedyness must be a float between 0 and 1.")
            return False

    if not os.path.exists(snakefile):
        logger.error("Error: Snakefile \"{}\" not present.".format(snakefile))
        return False
    snakefile = os.path.abspath(snakefile)

    if cluster and (drmaa is not None):
        raise ValueError("cluster and drmaa args are mutually exclusive")

    overwrite_config = dict()
    if configfile:
        overwrite_config.update(load_configfile(configfile))
    if config:
        overwrite_config.update(config)

    if workdir:
        olddir = os.getcwd()
        if not os.path.exists(workdir):
            logger.info("Creating specified working directory {}.".format(workdir))
            os.makedirs(workdir)
        workdir = os.path.abspath(workdir)
        os.chdir(workdir)
    workflow = Workflow(
        snakefile=snakefile, snakemakepath=snakemakepath,
        jobscript=jobscript, overwrite_shellcmd=overwrite_shellcmd,
        overwrite_config=overwrite_config, overwrite_workdir=workdir,
        overwrite_configfile=configfile,
        config_args=config_args
    )

    if standalone:
        try:
            # set the process group
            os.setpgrp()
        except:
            # ignore: if it does not work we can still work without it
            pass

    success = True
    try:
        workflow.include(
            snakefile,
            overwrite_first_rule=True,
            print_compilation=print_compilation
        )
        workflow.check()

        if not print_compilation:
            if listrules:
                workflow.list_rules()
            elif list_target_rules:
                workflow.list_rules(only_targets=True)
            elif list_resources:
                workflow.list_resources()
            else:
                    # if not printdag and not printrulegraph:
                    # handle subworkflows
                    subsnakemake = partial(
                        snakemake,
                        cores=cores,
                        nodes=nodes,
                        resources=resources,
                        dryrun=dryrun,
                        touch=touch,
                        printreason=printreason,
                        printshellcmds=printshellcmds,
                        nocolor=nocolor,
                        quiet=quiet,
                        keepgoing=keepgoing,
                        cluster=cluster,
                        drmaa=drmaa,
                        jobname=jobname,
                        immediate_submit=immediate_submit,
                        standalone=standalone,
                        ignore_ambiguity=ignore_ambiguity,
                        snakemakepath=snakemakepath,
                        lock=lock,
                        unlock=unlock,
                        cleanup_metadata=cleanup_metadata,
                        force_incomplete=force_incomplete,
                        ignore_incomplete=ignore_incomplete,
                        latency_wait=latency_wait,
                        benchmark_repeats=benchmark_repeats,
                        debug=debug,
                        notemp=notemp,
                        nodeps=nodeps,
                        jobscript=jobscript,
                        timestamp=timestamp,
                        greedyness=greedyness,
                        overwrite_shellcmd=overwrite_shellcmd,
                        keep_logger=True
                    )
                    success = workflow.execute(
                        targets=targets, dryrun=dryrun, touch=touch,
                        cores=cores, nodes=nodes, forcetargets=forcetargets,
                        forceall=forceall, forcerun=forcerun,
                        prioritytargets=prioritytargets, quiet=quiet,
                        keepgoing=keepgoing, printshellcmds=printshellcmds,
                        printreason=printreason, printrulegraph=printrulegraph,
                        printdag=printdag, cluster=cluster, jobname=jobname,
                        drmaa=drmaa,
                        printd3dag=printd3dag,
                        immediate_submit=immediate_submit,
                        ignore_ambiguity=ignore_ambiguity,
                        stats=stats,
                        force_incomplete=force_incomplete,
                        ignore_incomplete=ignore_incomplete,
                        list_version_changes=list_version_changes,
                        list_code_changes=list_code_changes,
                        list_input_changes=list_input_changes,
                        list_params_changes=list_params_changes,
                        summary=summary,
                        latency_wait=latency_wait,
                        benchmark_repeats=benchmark_repeats,
                        wait_for_files=wait_for_files,
                        detailed_summary=detailed_summary,
                        nolock=not lock,
                        unlock=unlock,
                        resources=resources,
                        notemp=notemp,
                        nodeps=nodeps,
                        keep_target_files=keep_target_files,
                        cleanup_metadata=cleanup_metadata,
                        subsnakemake=subsnakemake,
                        updated_files=updated_files,
                        allowed_rules=allowed_rules,
                        greedyness=greedyness
                    )

    # BrokenPipeError is not present in Python 3.2, so lets wait until everbody uses > 3.2
    #except BrokenPipeError:
        # ignore this exception and stop. It occurs if snakemake output is piped into less and less quits before reading the whole output.
        # in such a case, snakemake shall stop scheduling and quit with error 1
    #    success = False
    except (Exception, BaseException) as ex:
        print_exception(ex, workflow.linemaps)
        success = False
    if workdir:
        os.chdir(olddir)
    if workflow.persistence:
        workflow.persistence.unlock()
    if not keep_logger:
        logger.cleanup()
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


def parse_config(args):
    parsers = [int, float, eval, str]
    config = dict()
    if args.config is not None:
        valid = re.compile("[a-zA-Z_]\w*$")
        for entry in args.config:
            try:
                key, val = entry.split("=", 1)
            except ValueError:
                raise ValueError("Config entries have to be defined as name=value pairs.")
            if not valid.match(key):
                raise ValueError("Config entry must start with a valid identifier.")
            v = None
            for parser in parsers:
                try:
                    v = parser(val)
                    break
                except:
                    pass
            assert v is not None
            config[key] = v
    return config


def get_argument_parser():
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
        "--gui", nargs="?", const="8000", metavar="PORT", type=int,
        help="Serve an HTML based user interface to the given port "
        "(default: 8000). If possible, a browser window is opened.")
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
        "--config", nargs="*", metavar="KEY=VALUE",
        help=(
            "Set or overwrite values in the workflow config object. "
            "The workflow config object is accessible as variable config inside "
            "the workflow. Default values can be set by providing a JSON file "
            "(see Documentation)."))
    parser.add_argument(
        "--configfile", metavar="JSON_FILE",
        help=(
            "Specify or overwrite the config file of the workflow (see the docs). "
            "Values specified in JSON format are available in the global config "
            "dictionary inside the workflow."
        )
    )
    parser.add_argument(
        "--list", "-l", action="store_true",
        help="Show availiable rules in given Snakefile.")
    parser.add_argument(
        "--list-target-rules", "--lt", action="store_true",
        help="Show available target rules in given Snakefile.")
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
        "--rulegraph", action="store_true",
        help="Do not execute anything and print the dependency graph "
            "of rules in the dot language. This will be less "
            "crowded than above DAG of jobs, but also show less information. "
            "Note that each rule is displayed once, hence the displayed graph will be "
            "cyclic if a rule appears in several steps of the workflow. "
            "Use this if above option leads to a DAG that is too large. "
            "Recommended use on Unix systems: snakemake --ruledag | dot | display")
    parser.add_argument(
        "--d3dag", action="store_true",
        help="Print the DAG in D3.js compatible JSON format.")
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
        "--detailed-summary", "-D", action="store_true",
        help="Print a summary of all files created by the workflow. The "
        "has the following columns: filename, modification time, "
        "rule version, input file(s), shell command, status, plan.\n"
        "Thereby rule version contains the version"
        "the file was created with (see the version keyword of rules), and "
        "status denotes whether the file is missing, its input files are "
        "newer or if version or implementation of the rule changed since "
        "file creation. The input file and shell command columns are self"
        "explanatory. Finally the last column denotes whether the file "
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
            "files for a particular job are present.\n"
            "The submit command can be decorated to make it aware of certain job properties (input, output, params, wildcards, log, threads and dependencies (see the argument below)), e.g.:\n"
            "$ snakemake --cluster 'qsub -pe threaded {threads}'."))
    parser.add_argument(
        "--drmaa", nargs="?", const="", metavar="ARGS",
        help="Execute snakemake on a cluster accessed via DRMAA, "
            "Snakemake compiles jobs into scripts that are "
            "submitted to the cluster with the given command, once all input "
            "files for a particular job are present. ARGS can be used to "
            "specify options of the underlying cluster system, "
            "thereby using the job properties input, output, params, wildcards, log, "
            "threads and dependencies, e.g.: "
            "--drmaa ' -pe threaded {threads}'. Note that ARGS must be given in quotes and "
            "with a leading whitespace.")
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
        "--jobname", "--jn", default="snakejob.{rulename}.{jobid}.sh", metavar="NAME",
        help="Provide a custom name for the jobscript that is submitted to the cluster (see --cluster)."
        "NAME is \"snakejob.{rulename}.{jobid}.sh\" per default. The wildcard {jobid} has to be present in the name.")
    parser.add_argument(
        "--reason", "-r", action = "store_true",
        help="Print the reason for each executed rule.")
    parser.add_argument(
        "--stats", metavar="FILE",
        help="Write stats about Snakefile execution in JSON format to the given file.")
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
        help="List all output files that have been created with "
        "a different version (as determined by the version keyword).")
    parser.add_argument(
        "--list-code-changes", "--lc", action="store_true",
        help="List all output files for which the rule body (run or shell) have changed "
        "in the Snakefile.")
    parser.add_argument(
        "--list-input-changes", "--li", action="store_true",
        help="List all output files for which the defined input files have changed "
        "in the Snakefile (e.g. new input files were added in the rule definition or files were renamed). For listing input file modification in the filesystem, use --summary.")
    parser.add_argument(
        "--list-params-changes", "--lp", action="store_true",
        help="List all output files for which the defined params have changed "
        "in the Snakefile.")
    parser.add_argument(
        "--latency-wait", "--output-wait", "-w", type=int, default=5, metavar="SECONDS",
        help="Wait given seconds if an output file of a job is not present after "
        "the job finished. This helps if your filesystem "
        "suffers from latency (default 5).")
    parser.add_argument(
        "--wait-for-files", nargs="*", metavar="FILE", help="Wait --latency-wait seconds for these "
        "files to be present before executing the workflow. "
        "This option is used internally to handle filesystem latency in cluster "
        "environments.")
    parser.add_argument(
        "--benchmark-repeats", type=int, default=1, metavar="N",
        help="Repeat a job N times if marked for benchmarking (default 1)."
    )
    parser.add_argument(
        "--notemp", "--nt", action="store_true",
        help="Ignore temp() declarations. This is useful when running only "
        "a part of the workflow, since temp() would lead to deletion of "
        "probably needed files by other parts of the workflow."
        )
    parser.add_argument(
        "--keep-target-files", action="store_true",
        help="Do not adjust the paths of given target files relative to the working directory.")
    parser.add_argument(
        "--allowed-rules", nargs="+",
        help="Only use given rules. If omitted, all rules in Snakefile are used.")
    parser.add_argument(
        '--timestamp', '-T', action='store_true',
        help='Add a timestamp to all logging output')
    parser.add_argument(
        "--greedyness", type=float, default=None, help="Set the greedyness of scheduling. This value between 0 and 1 determines how careful jobs are selected for execution. The default value (1.0) provides the best speed and still acceptable scheduling quality.")
    parser.add_argument(
        "--print-compilation", action="store_true",
        help="Print the python representation of the workflow.")
    parser.add_argument(
        "--overwrite-shellcmd",
        help="Provide a shell command that shall be executed instead of those given in the workflow. "
        "This is for debugging purposes only.")
    parser.add_argument(
        "--debug", action="store_true", help="Print debugging output.")
    parser.add_argument(
        "--profile", metavar="FILE", help="Profile Snakemake and write the output to FILE. This requires yappi to be installed.")
    parser.add_argument(
        "--bash-completion", action="store_true", help="Output code to register bash completion for snakemake. Put the following in your .bashrc (including the accents): `snakemake --bash-completion` or issue it in an open terminal session.")
    parser.add_argument(
        "--version", "-v", action="version", version=__version__)
    return parser


def main():
    parser = get_argument_parser()
    args = parser.parse_args()

    if args.bash_completion:
        print("complete -C snakemake-bash-completion snakemake")
        sys.exit(0)

    snakemakepath = sys.argv[0]

    try:
        resources = parse_resources(args)
        config = parse_config(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        print("", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.cluster and args.drmaa:
        print("--cluster and --drmaa are mutually exclusive", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.profile:
        import yappi
        yappi.start()

    _snakemake = partial(snakemake, args.snakefile, snakemakepath=snakemakepath)

    if args.gui is not None:
        try:
            import snakemake.gui as gui
        except ImportError:
            print(
                "Error: GUI needs Flask to be installed. Install "
                "with easy_install or contact your administrator.",
                file=sys.stderr)
            sys.exit(1)

        _logging.getLogger("werkzeug").setLevel(_logging.ERROR)
        gui.register(_snakemake, args)
        url = "http://127.0.0.1:{}".format(args.gui)
        print("Listening on {}.".format(url), file=sys.stderr)
        def open_browser():
            try:
                webbrowser.open(url)
            except:
                pass
        print("Open this address in your browser to access the GUI.", file=sys.stderr)
        threading.Timer(0.5, open_browser).start()
        success = True
        try:
            gui.app.run(debug=False, threaded=True, port=args.gui)
        except (KeyboardInterrupt, SystemExit):
            # silently close
            pass
    else:
        success = _snakemake(
            listrules=args.list,
            list_target_rules=args.list_target_rules,
            cores=args.cores,
            nodes=args.cores,
            resources=resources,
            config=config,
            configfile=args.configfile,
            config_args=args.config,
            workdir=args.directory,
            targets=args.target,
            dryrun=args.dryrun,
            printshellcmds=args.printshellcmds,
            printreason=args.reason,
            printdag=args.dag,
            printrulegraph=args.rulegraph,
            printd3dag=args.d3dag,
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
            drmaa=args.drmaa,
            jobname=args.jobname,
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
            detailed_summary=args.detailed_summary,
            print_compilation=args.print_compilation,
            debug=args.debug,
            jobscript=args.jobscript,
            notemp=args.notemp,
            timestamp=args.timestamp,
            greedyness=args.greedyness,
            overwrite_shellcmd=args.overwrite_shellcmd,
            latency_wait=args.latency_wait,
            benchmark_repeats=args.benchmark_repeats,
            wait_for_files=args.wait_for_files,
            keep_target_files=args.keep_target_files,
            allowed_rules=args.allowed_rules)

    if args.profile:
        with open(args.profile, "w") as out:
            profile = yappi.get_func_stats()
            profile.sort("totaltime")
            profile.print_all(out=out)

    sys.exit(0 if success else 1)


def bash_completion(snakefile="Snakefile"):
    if not len(sys.argv) >= 2:
        print("Calculate bash completion for snakemake. This tool shall not be invoked by hand.")
        sys.exit(1)

    prefix = sys.argv[2]

    if prefix.startswith("-"):
        opts = [
            action.option_strings[0] for action in get_argument_parser()._actions
            if action.option_strings and action.option_strings[0].startswith(prefix)]
        print(*opts, sep="\n")
    else:
        files = glob.glob("{}*".format(prefix))
        if files:
            print(*files, sep="\n")
        elif os.path.exists(snakefile):
            workflow = Workflow(snakefile=snakefile, snakemakepath="snakemake")
            workflow.include(snakefile)

            workflow_files = sorted(set(
                file for file in workflow.concrete_files
                if file.startswith(prefix)))
            if workflow_files:
                print(*workflow_files, sep="\n")
            
            rules = [
                rule.name for rule in workflow.rules
                if rule.name.startswith(prefix)]
            if rules:
                print(*rules, sep="\n")
    sys.exit(0)
