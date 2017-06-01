__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import subprocess
import glob
import argparse
from argparse import ArgumentError
import logging as _logging
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
from snakemake.utils import update_config, available_cpu_count
from snakemake.common import Mode

def snakemake(snakefile,
              listrules=False,
              list_target_rules=False,
              cores=1,
              nodes=1,
              local_cores=1,
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
              until=[],
              omit_from=[],
              prioritytargets=[],
              stats=None,
              printreason=False,
              printshellcmds=False,
              debug_dag=False,
              printdag=False,
              printrulegraph=False,
              printd3dag=False,
              nocolor=False,
              quiet=False,
              keepgoing=False,
              cluster=None,
              cluster_config=None,
              cluster_sync=None,
              drmaa=None,
              drmaa_log_dir=None,
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
              archive=None,
              detailed_summary=False,
              latency_wait=3,
              benchmark_repeats=1,
              wait_for_files=None,
              print_compilation=False,
              debug=False,
              notemp=False,
              keep_remote_local=False,
              nodeps=False,
              keep_target_files=False,
              keep_shadow=False,
              allowed_rules=None,
              jobscript=None,
              timestamp=False,
              greediness=None,
              no_hooks=False,
              overwrite_shellcmd=None,
              updated_files=None,
              log_handler=None,
              keep_logger=False,
              max_jobs_per_second=None,
              restart_times=0,
              verbose=False,
              force_use_threads=False,
              use_conda=False,
              conda_prefix=None,
              mode=Mode.default,
              wrapper_prefix=None):
    """Run snakemake on a given snakefile.

    This function provides access to the whole snakemake functionality. It is not thread-safe.

    Args:
        snakefile (str):            the path to the snakefile
        listrules (bool):           list rules (default False)
        list_target_rules (bool):   list target rules (default False)
        cores (int):                the number of provided cores (ignored when using cluster support) (default 1)
        nodes (int):                the number of provided cluster nodes (ignored without cluster support) (default 1)
        local_cores (int):                the number of provided local cores if in cluster mode (ignored without cluster support) (default 1)
        resources (dict):           provided resources, a dictionary assigning integers to resource names, e.g. {gpu=1, io=5} (default {})
        config (dict):              override values for workflow config
        workdir (str):              path to working directory (default None)
        targets (list):             list of targets, e.g. rule or file names (default None)
        dryrun (bool):              only dry-run the workflow (default False)
        touch (bool):               only touch all output files if present (default False)
        forcetargets (bool):        force given targets to be re-created (default False)
        forceall (bool):            force all output files to be re-created (default False)
        forcerun (list):            list of files and rules that shall be re-created/re-executed (default [])
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
        cluster_config (str,list):  configuration file for cluster options, or list thereof (default None)
        cluster_sync (str):         blocking cluster submission command (like SGE 'qsub -sync y')  (default None)
        drmaa (str):                if not None use DRMAA for cluster support, str specifies native args passed to the cluster when submitting a job
        drmaa_log_dir (str):        the path to stdout and stderr output of DRMAA jobs (default None)
        jobname (str):              naming scheme for cluster job scripts (default "snakejob.{rulename}.{jobid}.sh")
        immediate_submit (bool):    immediately submit all cluster jobs, regardless of dependencies (default False)
        standalone (bool):          kill all processes very rudely in case of failure (do not use this if you use this API) (default False) (deprecated)
        ignore_ambiguity (bool):    ignore ambiguous rules and always take the first possible one (default False)
        snakemakepath (str):        Deprecated parameter whose value is ignored. Do not use.
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
        archive (str):              archive workflow into the given tarball
        latency_wait (int):         how many seconds to wait for an output file to appear after the execution of a job, e.g. to handle filesystem latency (default 3)
        benchmark_repeats (int):    number of repeated runs of a job if declared for benchmarking (default 1)
        wait_for_files (list):      wait for given files to be present before executing the workflow
        list_resources (bool):      list resources used in the workflow (default False)
        summary (bool):             list summary of all output files and their status (default False). If no option  is specified a basic summary will be ouput. If 'detailed' is added as an option e.g --summary detailed, extra info about the input and shell commands will be included
        detailed_summary (bool):    list summary of all input and output files and their status (default False)
        print_compilation (bool):   print the compilation of the snakefile (default False)
        debug (bool):               allow to use the debugger within rules
        notemp (bool):              ignore temp file flags, e.g. do not delete output files marked as temp after use (default False)
        keep_remote_local (bool):   keep local copies of remote files (default False)
        nodeps (bool):              ignore dependencies (default False)
        keep_target_files (bool):   Do not adjust the paths of given target files relative to the working directory.
        keep_shadow (bool):         Do not delete the shadow directory on snakemake startup.
        allowed_rules (set):        Restrict allowed rules to the given set. If None or empty, all rules are used.
        jobscript (str):            path to a custom shell script template for cluster jobs (default None)
        timestamp (bool):           print time stamps in front of any output (default False)
        greediness (float):         set the greediness of scheduling. This value between 0 and 1 determines how careful jobs are selected for execution. The default value (0.5 if prioritytargets are used, 1.0 else) provides the best speed and still acceptable scheduling quality.
        overwrite_shellcmd (str):   a shell command that shall be executed instead of those given in the workflow. This is for debugging purposes only.
        updated_files(list):        a list that will be filled with the files that are updated or created during the workflow execution
        verbose (bool):             show additional debug output (default False)
        max_jobs_per_second (int):  maximal number of cluster/drmaa jobs per second, None to impose no limit (default None)
        restart_times (int):        number of times to restart failing jobs (default 1)
        force_use_threads:          whether to force use of threads over processes. helpful if shared memory is full or unavailable (default False)
        use_conda (bool):           create conda environments for each job (defined with conda directive of rules)
        conda_prefix (str):            the directories in which conda environments will be created (default None)
        mode (snakemake.common.Mode): Execution mode
        wrapper_prefix (str):       Prefix for wrapper script URLs (default None)
        log_handler (function):     redirect snakemake output to this custom log handler, a function that takes a log message dictionary (see below) as its only argument (default None). The log message dictionary for the log handler has to following entries:

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
    assert not immediate_submit or (immediate_submit and notemp), "immediate_submit has to be combined with notemp (it does not support temp file handling)"

    if updated_files is None:
        updated_files = list()

    if cluster or cluster_sync or drmaa:
        cores = sys.maxsize
    else:
        nodes = sys.maxsize

    if isinstance(cluster_config, str):
        # Loading configuration from one file is still supported for
        # backward compatibility
        cluster_config = [cluster_config]
    if cluster_config:
        # Load all configuration files
        configs = [load_configfile(f) for f in cluster_config]
        # Merge in the order as specified, overriding earlier values with
        # later ones
        cluster_config = configs[0]
        for other in configs[1:]:
            update_config(cluster_config, other)
    else:
        cluster_config = dict()

    # force thread use for any kind of cluster
    use_threads = force_use_threads or (os.name != "posix") or cluster or cluster_sync or drmaa
    if not keep_logger:
        stdout = (
            (dryrun and not (printdag or printd3dag or printrulegraph)) or
            listrules or list_target_rules or list_resources
        )
        setup_logger(handler=log_handler,
                     quiet=quiet,
                     printreason=printreason,
                     printshellcmds=printshellcmds,
                     debug_dag=debug_dag,
                     nocolor=nocolor,
                     stdout=stdout,
                     debug=verbose,
                     timestamp=timestamp,
                     use_threads=use_threads,
                     mode=mode)

    if greediness is None:
        greediness = 0.5 if prioritytargets else 1.0
    else:
        if not (0 <= greediness <= 1.0):
            logger.error("Error: greediness must be a float between 0 and 1.")
            return False

    if not os.path.exists(snakefile):
        logger.error("Error: Snakefile \"{}\" not present.".format(snakefile))
        return False
    snakefile = os.path.abspath(snakefile)

    cluster_mode = (cluster is not None) + (cluster_sync is not
                                            None) + (drmaa is not None)
    if cluster_mode > 1:
        logger.error("Error: cluster and drmaa args are mutually exclusive")
        return False
    if debug and (cores > 1 or cluster_mode):
        logger.error(
            "Error: debug mode cannot be used with more than one core or cluster execution.")
        return False

    overwrite_config = dict()
    if configfile:
        overwrite_config.update(load_configfile(configfile))
        configfile = os.path.abspath(configfile)
    if config:
        overwrite_config.update(config)

    if workdir:
        olddir = os.getcwd()
        if not os.path.exists(workdir):
            logger.info(
                "Creating specified working directory {}.".format(workdir))
            os.makedirs(workdir)
        workdir = os.path.abspath(workdir)
        os.chdir(workdir)

    workflow = Workflow(snakefile=snakefile,
                        jobscript=jobscript,
                        overwrite_shellcmd=overwrite_shellcmd,
                        overwrite_config=overwrite_config,
                        overwrite_workdir=workdir,
                        overwrite_configfile=configfile,
                        overwrite_clusterconfig=cluster_config,
                        config_args=config_args,
                        debug=debug,
                        use_conda=use_conda,
                        conda_prefix=conda_prefix,
                        mode=mode,
                        wrapper_prefix=wrapper_prefix,
                        printshellcmds=printshellcmds,
                        restart_times=restart_times)
    success = True
    try:
        workflow.include(snakefile,
                         overwrite_first_rule=True,
                         print_compilation=print_compilation)
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
                subsnakemake = partial(snakemake,
                                       cores=cores,
                                       nodes=nodes,
                                       local_cores=local_cores,
                                       resources=resources,
                                       dryrun=dryrun,
                                       touch=touch,
                                       printreason=printreason,
                                       printshellcmds=printshellcmds,
                                       debug_dag=debug_dag,
                                       nocolor=nocolor,
                                       quiet=quiet,
                                       keepgoing=keepgoing,
                                       cluster=cluster,
                                       cluster_sync=cluster_sync,
                                       drmaa=drmaa,
                                       drmaa_log_dir=drmaa_log_dir,
                                       jobname=jobname,
                                       immediate_submit=immediate_submit,
                                       standalone=standalone,
                                       ignore_ambiguity=ignore_ambiguity,
                                       lock=lock,
                                       unlock=unlock,
                                       cleanup_metadata=cleanup_metadata,
                                       force_incomplete=force_incomplete,
                                       ignore_incomplete=ignore_incomplete,
                                       latency_wait=latency_wait,
                                       benchmark_repeats=benchmark_repeats,
                                       verbose=verbose,
                                       notemp=notemp,
                                       keep_remote_local=keep_remote_local,
                                       nodeps=nodeps,
                                       jobscript=jobscript,
                                       timestamp=timestamp,
                                       greediness=greediness,
                                       no_hooks=no_hooks,
                                       overwrite_shellcmd=overwrite_shellcmd,
                                       config=config,
                                       config_args=config_args,
                                       cluster_config=cluster_config,
                                       keep_logger=True,
                                       keep_shadow=True,
                                       force_use_threads=use_threads,
                                       use_conda=use_conda,
                                       conda_prefix=conda_prefix)
                success = workflow.execute(
                    targets=targets,
                    dryrun=dryrun,
                    touch=touch,
                    cores=cores,
                    nodes=nodes,
                    local_cores=local_cores,
                    forcetargets=forcetargets,
                    forceall=forceall,
                    forcerun=forcerun,
                    prioritytargets=prioritytargets,
                    until=until,
                    omit_from=omit_from,
                    quiet=quiet,
                    keepgoing=keepgoing,
                    printshellcmds=printshellcmds,
                    printreason=printreason,
                    printrulegraph=printrulegraph,
                    printdag=printdag,
                    cluster=cluster,
                    cluster_sync=cluster_sync,
                    jobname=jobname,
                    drmaa=drmaa,
                    drmaa_log_dir=drmaa_log_dir,
                    max_jobs_per_second=max_jobs_per_second,
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
                    archive=archive,
                    latency_wait=latency_wait,
                    benchmark_repeats=benchmark_repeats,
                    wait_for_files=wait_for_files,
                    detailed_summary=detailed_summary,
                    nolock=not lock,
                    unlock=unlock,
                    resources=resources,
                    notemp=notemp,
                    keep_remote_local=keep_remote_local,
                    nodeps=nodeps,
                    keep_target_files=keep_target_files,
                    keep_shadow=keep_shadow,
                    cleanup_metadata=cleanup_metadata,
                    subsnakemake=subsnakemake,
                    updated_files=updated_files,
                    allowed_rules=allowed_rules,
                    greediness=greediness,
                    no_hooks=no_hooks,
                    force_use_threads=use_threads)

    except BrokenPipeError:
        # ignore this exception and stop. It occurs if snakemake output is piped into less and less quits before reading the whole output.
        # in such a case, snakemake shall stop scheduling and quit with error 1
        success = False
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
    """Parse resources from args."""
    resources = dict()
    if args.resources is not None:
        valid = re.compile("[a-zA-Z_]\w*$")
        for res in args.resources:
            try:
                res, val = res.split("=")
            except ValueError:
                raise ValueError(
                    "Resources have to be defined as name=value pairs.")
            if not valid.match(res):
                raise ValueError(
                    "Resource definition must start with a valid identifier.")
            try:
                val = int(val)
            except ValueError:
                raise ValueError(
                    "Resource definiton must contain an integer after the identifier.")
            if res == "_cores":
                raise ValueError(
                    "Resource _cores is already defined internally. Use a different name.")
            resources[res] = val
    return resources


def parse_config(args):
    """Parse config from args."""
    parsers = [int, float, eval, str]
    config = dict()
    if args.config is not None:
        valid = re.compile("[a-zA-Z_]\w*$")
        for entry in args.config:
            try:
                key, val = entry.split("=", 1)
            except ValueError:
                raise ValueError(
                    "Config entries have to be defined as name=value pairs.")
            if not valid.match(key):
                raise ValueError(
                    "Config entry must start with a valid identifier.")
            v = None
            for parser in parsers:
                try:
                    v = parser(val)
                    # avoid accidental interpretation as function
                    if not callable(v):
                        break
                except:
                    pass
            assert v is not None
            config[key] = v
    return config


def get_argument_parser():
    """Generate and return argument parser."""
    parser = argparse.ArgumentParser(
        description="Snakemake is a Python based language and execution "
        "environment for GNU Make-like workflows.")

    parser.add_argument("target",
                        nargs="*",
                        default=None,
                        help="Targets to build. May be rules or files.")
    parser.add_argument("--snakefile", "-s",
                        metavar="FILE",
                        default="Snakefile",
                        help="The workflow definition in a snakefile.")
    parser.add_argument(
        "--gui",
        nargs="?",
        const="8000",
        metavar="PORT",
        type=int,
        help="Serve an HTML based user interface to the given port "
        "(default: 8000). If possible, a browser window is opened.")
    parser.add_argument(
        "--cores", "--jobs", "-j",
        action="store",
        const=available_cpu_count(),
        nargs="?",
        metavar="N",
        type=int,
        help=("Use at most N cores in parallel (default: 1). "
              "If N is omitted, the limit is set to the number of "
              "available cores."))
    parser.add_argument(
        "--local-cores",
        action="store",
        default=available_cpu_count(),
        metavar="N",
        type=int,
        help=
        ("In cluster mode, use at most N cores of the host machine in parallel "
         " (default: number of CPU cores of the host). The cores are used to execute "
         "local rules. This option is ignored when not in cluster mode."))
    parser.add_argument(
        "--resources", "--res",
        nargs="*",
        metavar="NAME=INT",
        help=("Define additional resources that shall constrain the scheduling "
              "analogously to threads (see above). A resource is defined as "
              "a name and an integer value. E.g. --resources gpu=1. Rules can "
              "use resources by defining the resource keyword, e.g. "
              "resources: gpu=1. If now two rules require 1 of the resource "
              "'gpu' they won't be run in parallel by the scheduler."))
    parser.add_argument(
        "--config", "-C",
        nargs="*",
        metavar="KEY=VALUE",
        help=
        ("Set or overwrite values in the workflow config object. "
         "The workflow config object is accessible as variable config inside "
         "the workflow. Default values can be set by providing a JSON file "
         "(see Documentation)."))
    parser.add_argument(
        "--configfile",
        metavar="FILE",
        help=
        ("Specify or overwrite the config file of the workflow (see the docs). "
         "Values specified in JSON or YAML format are available in the global config "
         "dictionary inside the workflow."))
    parser.add_argument("--list", "-l",
                        action="store_true",
                        help="Show availiable rules in given Snakefile.")
    parser.add_argument("--list-target-rules", "--lt",
                        action="store_true",
                        help="Show available target rules in given Snakefile.")
    parser.add_argument("--directory", "-d",
                        metavar="DIR",
                        action="store",
                        help=("Specify working directory (relative paths in "
                              "the snakefile will use this as their origin)."))
    parser.add_argument("--dryrun", "-n",
                        action="store_true",
                        help="Do not execute anything.")
    parser.add_argument(
        "--printshellcmds", "-p",
        action="store_true",
        help="Print out the shell commands that will be executed.")
    parser.add_argument(
        "--debug-dag",
        action="store_true",
        help="Print candidate and selected jobs (including their wildcards) while "
        "inferring DAG. This can help to debug unexpected DAG topology or errors.")
    parser.add_argument(
        "--dag",
        action="store_true",
        help="Do not execute anything and print the directed "
        "acyclic graph of jobs in the dot language. Recommended "
        "use on Unix systems: snakemake --dag | dot | display")
    parser.add_argument(
        "--force-use-threads",
        dest="force_use_threads",
        action="store_true",
        help="Force threads rather than processes. Helpful if shared memory (/dev/shm) is full or unavailable.")
    parser.add_argument(
        "--rulegraph",
        action="store_true",
        help="Do not execute anything and print the dependency graph "
        "of rules in the dot language. This will be less "
        "crowded than above DAG of jobs, but also show less information. "
        "Note that each rule is displayed once, hence the displayed graph will be "
        "cyclic if a rule appears in several steps of the workflow. "
        "Use this if above option leads to a DAG that is too large. "
        "Recommended use on Unix systems: snakemake --rulegraph | dot | display")
    parser.add_argument("--d3dag",
                        action="store_true",
                        help="Print the DAG in D3.js compatible JSON format.")
    parser.add_argument(
        "--summary", "-S",
        action="store_true",
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
        "--detailed-summary", "-D",
        action="store_true",
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
        "--archive",
        metavar="FILE",
        help="Archive the workflow into the given tar archive FILE. The archive "
        "will be created such that the workflow can be re-executed on a vanilla "
        "system. The function needs conda and git to be installed. "
        "It will archive every file that is under git version control. "
        "Note that it is best practice to have the Snakefile, config files, and "
        "scripts under version control. Hence, they will be included in the archive. "
        "Further, it will add input files that are not generated by "
        "by the workflow itself and conda environments. Note that symlinks are "
        "dereferenced. Supported "
        "formats are .tar, .tar.gz, .tar.bz2 and .tar.xz."
    )
    parser.add_argument(
        "--touch", "-t",
        action="store_true",
        help=("Touch output files (mark them up to date without really "
              "changing them) instead of running their commands. This is "
              "used to pretend that the rules were executed, in order to "
              "fool future invocations of snakemake. Fails if a file does "
              "not yet exist."))
    parser.add_argument("--keep-going", "-k",
                        action="store_true",
                        help="Go on with independent jobs if a job fails.")
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help=("Force the execution of the selected target or the first rule "
              "regardless of already created output."))
    parser.add_argument(
        "--forceall", "-F",
        action="store_true",
        help=("Force the execution of the selected (or the first) rule and "
              "all rules it is dependent on regardless of already created "
              "output."))
    parser.add_argument(
        "--forcerun", "-R",
        nargs="*",
        metavar="TARGET",
        help=("Force the re-execution or creation of the given rules or files."
              " Use this option if you changed a rule and want to have all its "
              "output in your workflow updated."))
    parser.add_argument(
        "--prioritize", "-P",
        nargs="+",
        metavar="TARGET",
        help=("Tell the scheduler to assign creation of given targets "
              "(and all their dependencies) highest priority. (EXPERIMENTAL)"))
    parser.add_argument(
        "--until", "-U",
        nargs="+",
        metavar="TARGET",
        help=("Runs the pipeline until it reaches the specified rules or "
              "files. Only runs jobs that are dependencies of the specified "
              "rule or files, does not run sibling DAGs. "))
    parser.add_argument(
        "--omit-from", "-O",
        nargs="+",
        metavar="TARGET",
        help=("Prevent the execution or creation of the given rules or files "
              "as well as any rules or files that are downstream of these targets "
              "in the DAG. Also runs jobs in sibling DAGs that are independent of the "
              "rules or files specified here."))
    parser.add_argument(
        "--allow-ambiguity", "-a",
        action="store_true",
        help=("Don't check for ambiguous rules and simply use the first if "
              "several can produce the same file. This allows the user to "
              "prioritize rules by their order in the snakefile."))
    # TODO extend below description to explain the wildcards that can be used

    cluster_group = parser.add_mutually_exclusive_group()
    cluster_group.add_argument(
        "--cluster", "-c",
        metavar="CMD",
        help=
        ("Execute snakemake rules with the given submit command, "
         "e.g. qsub. Snakemake compiles jobs into scripts that are "
         "submitted to the cluster with the given command, once all input "
         "files for a particular job are present.\n"
         "The submit command can be decorated to make it aware of certain "
         "job properties (input, output, params, wildcards, log, threads "
         "and dependencies (see the argument below)), e.g.:\n"
         "$ snakemake --cluster 'qsub -pe threaded {threads}'.")),
    cluster_group.add_argument(
        "--cluster-sync",
        metavar="CMD",
        help=
        ("cluster submission command will block, returning the remote exit"
         "status upon remote termination (for example, this should be used"
         "if the cluster command is 'qsub -sync y' (SGE)")),
    cluster_group.add_argument(
        "--drmaa",
        nargs="?",
        const="",
        metavar="ARGS",
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
        "--drmaa-log-dir",
        metavar="DIR",
        help="Specify a directory in which stdout and stderr files of DRMAA"
        " jobs will be written. The value may be given as a relative path,"
        " in which case Snakemake will use the current invocation directory"
        " as the origin. If given, this will override any given '-o' and/or"
        " '-e' native specification. If not given, all DRMAA stdout and"
        " stderr files are written to the current working directory.")

    parser.add_argument(
        "--cluster-config", "-u",
        metavar="FILE",
        default=[],
        action="append",
        help=
        ("A JSON or YAML file that defines the wildcards used in 'cluster'"
         "for specific rules, instead of having them specified in the Snakefile. "
         "For example, for rule 'job' you may define: "
         "{ 'job' : { 'time' : '24:00:00' } } to specify the time for rule 'job'. "
         "You can specify more than one file.  The configuration files are merged "
         "with later values overriding earlier ones.")),
    parser.add_argument(
        "--immediate-submit", "--is",
        action="store_true",
        help="Immediately submit all jobs to the cluster instead of waiting "
        "for present input files. This will fail, unless you make "
        "the cluster aware of job dependencies, e.g. via:\n"
        "$ snakemake --cluster 'sbatch --dependency {dependencies}.\n"
        "Assuming that your submit script (here sbatch) outputs the "
        "generated job id to the first stdout line, {dependencies} will "
        "be filled with space separated job ids this job depends on.")
    parser.add_argument(
        "--jobscript", "--js",
        metavar="SCRIPT",
        help="Provide a custom job script for submission to the cluster. "
        "The default script resides as 'jobscript.sh' in the "
        "installation directory.")
    parser.add_argument(
        "--jobname", "--jn",
        default="snakejob.{rulename}.{jobid}.sh",
        metavar="NAME",
        help="Provide a custom name for the jobscript that is submitted to the "
        "cluster (see --cluster). NAME is \"snakejob.{rulename}.{jobid}.sh\" "
        "per default. The wildcard {jobid} has to be present in the name.")
    parser.add_argument("--reason", "-r",
                        action="store_true",
                        help="Print the reason for each executed rule.")
    parser.add_argument(
        "--stats",
        metavar="FILE",
        help=
        "Write stats about Snakefile execution in JSON format to the given file.")
    parser.add_argument("--nocolor",
                        action="store_true",
                        help="Do not use a colored output.")
    parser.add_argument("--quiet", "-q",
                        action="store_true",
                        help="Do not output any progress or rule information.")
    parser.add_argument("--nolock",
                        action="store_true",
                        help="Do not lock the working directory")
    parser.add_argument("--unlock",
                        action="store_true",
                        help="Remove a lock on the working directory.")
    parser.add_argument(
        "--cleanup-metadata", "--cm",
        nargs="+",
        metavar="FILE",
        help="Cleanup the metadata "
        "of given files. That means that snakemake removes any tracked "
        "version info, and any marks that files are incomplete.")
    parser.add_argument(
        "--rerun-incomplete", "--ri",
        action="store_true",
        help="Re-run all "
        "jobs the output of which is recognized as incomplete.")
    parser.add_argument("--ignore-incomplete", "--ii",
                        action="store_true",
                        help="Do not check for incomplete output files.")
    parser.add_argument(
        "--list-version-changes", "--lv",
        action="store_true",
        help="List all output files that have been created with "
        "a different version (as determined by the version keyword).")
    parser.add_argument(
        "--list-code-changes", "--lc",
        action="store_true",
        help=
        "List all output files for which the rule body (run or shell) have "
        "changed in the Snakefile.")
    parser.add_argument(
        "--list-input-changes", "--li",
        action="store_true",
        help=
        "List all output files for which the defined input files have changed "
        "in the Snakefile (e.g. new input files were added in the rule "
        "definition or files were renamed). For listing input file "
        "modification in the filesystem, use --summary.")
    parser.add_argument(
        "--list-params-changes", "--lp",
        action="store_true",
        help="List all output files for which the defined params have changed "
        "in the Snakefile.")
    parser.add_argument(
        "--latency-wait", "--output-wait", "-w",
        type=int,
        default=5,
        metavar="SECONDS",
        help=
        "Wait given seconds if an output file of a job is not present after "
        "the job finished. This helps if your filesystem "
        "suffers from latency (default 5).")
    parser.add_argument(
        "--wait-for-files",
        nargs="*",
        metavar="FILE",
        help="Wait --latency-wait seconds for these "
        "files to be present before executing the workflow. "
        "This option is used internally to handle filesystem latency in cluster "
        "environments.")
    parser.add_argument(
        "--benchmark-repeats",
        type=int,
        default=1,
        metavar="N",
        help="Repeat a job N times if marked for benchmarking (default 1).")
    parser.add_argument(
        "--notemp", "--nt",
        action="store_true",
        help="Ignore temp() declarations. This is useful when running only "
        "a part of the workflow, since temp() would lead to deletion of "
        "probably needed files by other parts of the workflow.")
    parser.add_argument(
        "--keep-remote",
        action="store_true",
        help="Keep local copies of remote input files.")
    parser.add_argument(
        "--keep-target-files",
        action="store_true",
        help=
        "Do not adjust the paths of given target files relative to the working directory.")
    parser.add_argument(
        "--keep-shadow",
        action="store_true",
        help=
        "Do not delete the shadow directory on snakemake startup.")
    parser.add_argument(
        "--allowed-rules",
        nargs="+",
        help=
        "Only consider given rules. If omitted, all rules in Snakefile are "
        "used. Note that this is intended primarily for internal use and may "
        "lead to unexpected results otherwise.")
    parser.add_argument(
        "--max-jobs-per-second", default=None, type=float,
        help=
        "Maximal number of cluster/drmaa jobs per second, default is no limit")
    parser.add_argument(
        "--restart-times", default=0, type=int,
        help=
        "Number of times to restart failing jobs (defaults to 0).")
    parser.add_argument('--timestamp', '-T',
                        action='store_true',
                        help='Add a timestamp to all logging output')
    parser.add_argument(
        "--greediness",
        type=float,
        default=None,
        help="Set the greediness of scheduling. This value between 0 and 1 "
        "determines how careful jobs are selected for execution. The default "
        "value (1.0) provides the best speed and still acceptable scheduling "
        "quality.")
    parser.add_argument(
        "--no-hooks",
        action="store_true",
        help="Do not invoke onstart, onsuccess or onerror hooks after execution.")
    parser.add_argument(
        "--print-compilation",
        action="store_true",
        help="Print the python representation of the workflow.")
    parser.add_argument(
        "--overwrite-shellcmd",
        help="Provide a shell command that shall be executed instead of those "
        "given in the workflow. "
        "This is for debugging purposes only.")
    parser.add_argument("--verbose",
                        action="store_true",
                        help="Print debugging output.")
    parser.add_argument("--debug",
                        action="store_true",
                        help="Allow to debug rules with e.g. PDB. This flag "
                        "allows to set breakpoints in run blocks.")
    parser.add_argument(
        "--profile",
        metavar="FILE",
        help=
        "Profile Snakemake and write the output to FILE. This requires yappi "
        "to be installed.")
    parser.add_argument(
        "--mode",
        choices=[Mode.default, Mode.subprocess, Mode.cluster],
        default=Mode.default,
        type=int,
        help="Set execution mode of Snakemake (internal use only)."
    )
    parser.add_argument(
        "--bash-completion",
        action="store_true",
        help="Output code to register bash completion for snakemake. Put the "
        "following in your .bashrc (including the accents): "
        "`snakemake --bash-completion` or issue it in an open terminal "
        "session.")
    parser.add_argument(
        "--use-conda",
        action="store_true",
        help="If defined in the rule, create job specific conda environments. "
        "If this flag is not set, the conda directive is ignored.")
    parser.add_argument(
        "--conda-prefix",
        metavar="DIR",
        help="Specify a directory in which the 'conda' and 'conda-archive' "
        "directories are created. These are used to store conda environments "
        "and their archives, respectively. If not supplied, the value is set "
        "to the '.snakemake' directory relative to the invocation directory. "
        "If supplied, the `--use-conda` flag must also be set. The value may "
        "be given as a relative path, which will be extrapolated to the "
        "invocation directory, or as an absolute path.")
    parser.add_argument(
        "--wrapper-prefix",
        default="https://bitbucket.org/snakemake/snakemake-wrappers/raw/",
        help="Prefix for URL created from wrapper directive (default: "
        "https://bitbucket.org/snakemake/snakemake-wrappers/raw/). Set this to "
        "a different URL to use your fork or a local clone of the repository."
    )
    parser.add_argument("--version", "-v",
                        action="version",
                        version=__version__)
    return parser


def main(argv=None):
    """Main entry point."""
    parser = get_argument_parser()
    args = parser.parse_args(argv)

    if args.bash_completion:
        cmd = b"complete -o bashdefault -C snakemake-bash-completion snakemake"
        sys.stdout.buffer.write(cmd)
        sys.exit(0)

    try:
        resources = parse_resources(args)
        config = parse_config(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        print("", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if (args.cluster or args.cluster_sync or args.drmaa):
        if args.cores is None:
            if args.dryrun:
                args.cores = 1
            else:
                print(
                    "Error: you need to specify the maximum number of jobs to "
                    "be queued or executed at the same time with --jobs.",
                    file=sys.stderr)
                sys.exit(1)
    elif args.cores is None:
        args.cores = 1

    if args.drmaa_log_dir is not None:
        if not os.path.isabs(args.drmaa_log_dir):
            args.drmaa_log_dir = os.path.abspath(os.path.expanduser(args.drmaa_log_dir))

    if args.profile:
        import yappi
        yappi.start()

    if args.immediate_submit and not args.notemp:
        print(
            "Error: --immediate-submit has to be combined with --notemp, "
            "because temp file handling is not supported in this mode.",
            file=sys.stderr)
        sys.exit(1)

    if args.conda_prefix and not args.use_conda:
        print(
            "Error: --use-conda must be set if --conda-prefix is set.",
            file=sys.stderr)
        sys.exit(1)

    if args.gui is not None:
        try:
            import snakemake.gui as gui
        except ImportError:
            print("Error: GUI needs Flask to be installed. Install "
                  "with easy_install or contact your administrator.",
                  file=sys.stderr)
            sys.exit(1)

        _logging.getLogger("werkzeug").setLevel(_logging.ERROR)

        _snakemake = partial(snakemake, os.path.abspath(args.snakefile))
        gui.register(_snakemake, args)
        url = "http://127.0.0.1:{}".format(args.gui)
        print("Listening on {}.".format(url), file=sys.stderr)

        def open_browser():
            try:
                webbrowser.open(url)
            except:
                pass

        print("Open this address in your browser to access the GUI.",
              file=sys.stderr)
        threading.Timer(0.5, open_browser).start()
        success = True
        try:
            gui.app.run(debug=False, threaded=True, port=args.gui)
        except (KeyboardInterrupt, SystemExit):
            # silently close
            pass
    else:
        success = snakemake(args.snakefile,
                            listrules=args.list,
                            list_target_rules=args.list_target_rules,
                            cores=args.cores,
                            local_cores=args.local_cores,
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
                            debug_dag=args.debug_dag,
                            printdag=args.dag,
                            printrulegraph=args.rulegraph,
                            printd3dag=args.d3dag,
                            touch=args.touch,
                            forcetargets=args.force,
                            forceall=args.forceall,
                            forcerun=args.forcerun,
                            prioritytargets=args.prioritize,
                            until=args.until,
                            omit_from=args.omit_from,
                            stats=args.stats,
                            nocolor=args.nocolor,
                            quiet=args.quiet,
                            keepgoing=args.keep_going,
                            cluster=args.cluster,
                            cluster_config=args.cluster_config,
                            cluster_sync=args.cluster_sync,
                            drmaa=args.drmaa,
                            drmaa_log_dir=args.drmaa_log_dir,
                            jobname=args.jobname,
                            immediate_submit=args.immediate_submit,
                            standalone=True,
                            ignore_ambiguity=args.allow_ambiguity,
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
                            archive=args.archive,
                            print_compilation=args.print_compilation,
                            verbose=args.verbose,
                            debug=args.debug,
                            jobscript=args.jobscript,
                            notemp=args.notemp,
                            keep_remote_local=args.keep_remote,
                            timestamp=args.timestamp,
                            greediness=args.greediness,
                            no_hooks=args.no_hooks,
                            overwrite_shellcmd=args.overwrite_shellcmd,
                            latency_wait=args.latency_wait,
                            benchmark_repeats=args.benchmark_repeats,
                            wait_for_files=args.wait_for_files,
                            keep_target_files=args.keep_target_files,
                            keep_shadow=args.keep_shadow,
                            allowed_rules=args.allowed_rules,
                            max_jobs_per_second=args.max_jobs_per_second,
                            restart_times=args.restart_times,
                            force_use_threads=args.force_use_threads,
                            use_conda=args.use_conda,
                            conda_prefix=args.conda_prefix,
                            mode=args.mode,
                            wrapper_prefix=args.wrapper_prefix)

    if args.profile:
        with open(args.profile, "w") as out:
            profile = yappi.get_func_stats()
            profile.sort("totaltime")
            profile.print_all(out=out)

    sys.exit(0 if success else 1)


def bash_completion(snakefile="Snakefile"):
    """Entry point for bash completion."""
    if not len(sys.argv) >= 2:
        print(
            "Calculate bash completion for snakemake. This tool shall not be invoked by hand.")
        sys.exit(1)

    def print_candidates(candidates):
        if candidates:
            candidates = sorted(set(candidates))
            ## Use bytes for avoiding '^M' under Windows.
            sys.stdout.buffer.write(b'\n'.join(s.encode() for s in candidates))

    prefix = sys.argv[2]

    if prefix.startswith("-"):
        print_candidates(action.option_strings[0]
                         for action in get_argument_parser()._actions
                         if action.option_strings and
                         action.option_strings[0].startswith(prefix))
    else:
        files = glob.glob("{}*".format(prefix))
        if files:
            print_candidates(files)
        elif os.path.exists(snakefile):
            workflow = Workflow(snakefile=snakefile)
            workflow.include(snakefile)

            print_candidates([file
                              for file in workflow.concrete_files
                              if file.startswith(prefix)] +
                             [rule.name
                              for rule in workflow.rules
                              if rule.name.startswith(prefix)])
    sys.exit(0)
