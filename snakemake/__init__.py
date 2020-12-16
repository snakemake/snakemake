__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import subprocess
import glob
from argparse import ArgumentError, ArgumentDefaultsHelpFormatter
import logging as _logging
import re
import sys
import inspect
import threading
import webbrowser
from functools import partial
import importlib
import shutil
from importlib.machinery import SourceFileLoader

from snakemake.workflow import Workflow
from snakemake.dag import Batch
from snakemake.exceptions import print_exception, WorkflowError
from snakemake.logging import setup_logger, logger, SlackLogger
from snakemake.io import load_configfile
from snakemake.shell import shell
from snakemake.utils import update_config, available_cpu_count
from snakemake.common import Mode, __version__, MIN_PY_VERSION
from snakemake.resources import parse_resources, DefaultResources


SNAKEFILE_CHOICES = [
    "Snakefile",
    "snakefile",
    "workflow/Snakefile",
    "workflow/snakefile",
]


def snakemake(
    snakefile,
    batch=None,
    cache=None,
    report=None,
    report_stylesheet=None,
    lint=None,
    generate_unit_tests=None,
    listrules=False,
    list_target_rules=False,
    cores=1,
    nodes=1,
    local_cores=1,
    resources=dict(),
    overwrite_threads=None,
    overwrite_scatter=None,
    default_resources=None,
    config=dict(),
    configfiles=None,
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
    printfilegraph=False,
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
    conda_cleanup_envs=False,
    cleanup_shadow=False,
    cleanup_scripts=True,
    force_incomplete=False,
    ignore_incomplete=False,
    list_version_changes=False,
    list_code_changes=False,
    list_input_changes=False,
    list_params_changes=False,
    list_untracked=False,
    list_resources=False,
    summary=False,
    archive=None,
    delete_all_output=False,
    delete_temp_output=False,
    detailed_summary=False,
    latency_wait=3,
    wait_for_files=None,
    print_compilation=False,
    debug=False,
    notemp=False,
    keep_remote_local=False,
    nodeps=False,
    keep_target_files=False,
    allowed_rules=None,
    jobscript=None,
    greediness=None,
    no_hooks=False,
    overwrite_shellcmd=None,
    updated_files=None,
    log_handler=[],
    keep_logger=False,
    wms_monitor=None,
    max_jobs_per_second=None,
    max_status_checks_per_second=100,
    restart_times=0,
    attempt=1,
    verbose=False,
    force_use_threads=False,
    use_conda=False,
    use_singularity=False,
    use_env_modules=False,
    singularity_args="",
    conda_frontend="conda",
    conda_prefix=None,
    conda_cleanup_pkgs=None,
    list_conda_envs=False,
    singularity_prefix=None,
    shadow_prefix=None,
    scheduler="ilp",
    scheduler_ilp_solver=None,
    conda_create_envs_only=False,
    mode=Mode.default,
    wrapper_prefix=None,
    kubernetes=None,
    container_image=None,
    tibanna=False,
    tibanna_sfn=None,
    google_lifesciences=False,
    google_lifesciences_regions=None,
    google_lifesciences_location=None,
    google_lifesciences_cache=False,
    tes=None,
    preemption_default=None,
    preemptible_rules=None,
    precommand="",
    default_remote_provider=None,
    default_remote_prefix="",
    tibanna_config=False,
    assume_shared_fs=True,
    cluster_status=None,
    export_cwl=None,
    show_failed_logs=False,
    keep_incomplete=False,
    keep_metadata=True,
    messaging=None,
    edit_notebook=None,
    envvars=None,
    overwrite_groups=None,
    group_components=None,
    max_inventory_wait_time=20,
    execute_subworkflows=True,
):
    """Run snakemake on a given snakefile.

    This function provides access to the whole snakemake functionality. It is not thread-safe.

    Args:
        snakefile (str):            the path to the snakefile
        batch (Batch):              whether to compute only a partial DAG, defined by the given Batch object (default None)
        report (str):               create an HTML report for a previous run at the given path
        lint (str):                 print lints instead of executing (None, "plain" or "json", default None)
        listrules (bool):           list rules (default False)
        list_target_rules (bool):   list target rules (default False)
        cores (int):                the number of provided cores (ignored when using cluster support) (default 1)
        nodes (int):                the number of provided cluster nodes (ignored without cluster support) (default 1)
        local_cores (int):          the number of provided local cores if in cluster mode (ignored without cluster support) (default 1)
        resources (dict):           provided resources, a dictionary assigning integers to resource names, e.g. {gpu=1, io=5} (default {})
        default_resources (DefaultResources):   default values for resources not defined in rules (default None)
        config (dict):              override values for workflow config
        workdir (str):              path to working directory (default None)
        targets (list):             list of targets, e.g. rule or file names (default None)
        dryrun (bool):              only dry-run the workflow (default False)
        touch (bool):               only touch all output files if present (default False)
        forcetargets (bool):        force given targets to be re-created (default False)
        forceall (bool):            force all output files to be re-created (default False)
        forcerun (list):            list of files and rules that shall be re-created/re-executed (default [])
        execute_subworkflows (bool):   execute subworkflows if present (default True)
        prioritytargets (list):     list of targets that shall be run with maximum priority (default [])
        stats (str):                path to file that shall contain stats about the workflow execution (default None)
        printreason (bool):         print the reason for the execution of each job (default false)
        printshellcmds (bool):      print the shell command of each job (default False)
        printdag (bool):            print the dag in the graphviz dot language (default False)
        printrulegraph (bool):      print the graph of rules in the graphviz dot language (default False)
        printfilegraph (bool):      print the graph of rules with their input and output files in the graphviz dot language (default False)
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
        snakemakepath (str):        deprecated parameter whose value is ignored. Do not use.
        lock (bool):                lock the working directory when executing the workflow (default True)
        unlock (bool):              just unlock the working directory (default False)
        cleanup_metadata (list):    just cleanup metadata of given list of output files (default None)
        drop_metadata (bool):       drop metadata file tracking information after job finishes (--report and --list_x_changes information will be incomplete) (default False)
        conda_cleanup_envs (bool):  just cleanup unused conda environments (default False)
        cleanup_shadow (bool):      just cleanup old shadow directories (default False)
        cleanup_scripts (bool):     delete wrapper scripts used for execution (default True)
        force_incomplete (bool):    force the re-creation of incomplete files (default False)
        ignore_incomplete (bool):   ignore incomplete files (default False)
        list_version_changes (bool): list output files with changed rule version (default False)
        list_code_changes (bool):   list output files with changed rule code (default False)
        list_input_changes (bool):  list output files with changed input files (default False)
        list_params_changes (bool): list output files with changed params (default False)
        list_untracked (bool):      list files in the workdir that are not used in the workflow (default False)
        summary (bool):             list summary of all output files and their status (default False)
        archive (str):              archive workflow into the given tarball
        delete_all_output (bool)    remove all files generated by the workflow (default False)
        delete_temp_output (bool)   remove all temporary files generated by the workflow (default False)
        latency_wait (int):         how many seconds to wait for an output file to appear after the execution of a job, e.g. to handle filesystem latency (default 3)
        wait_for_files (list):      wait for given files to be present before executing the workflow
        list_resources (bool):      list resources used in the workflow (default False)
        summary (bool):             list summary of all output files and their status (default False). If no option  is specified a basic summary will be ouput. If 'detailed' is added as an option e.g --summary detailed, extra info about the input and shell commands will be included
        detailed_summary (bool):    list summary of all input and output files and their status (default False)
        print_compilation (bool):   print the compilation of the snakefile (default False)
        debug (bool):               allow to use the debugger within rules
        notemp (bool):              ignore temp file flags, e.g. do not delete output files marked as temp after use (default False)
        keep_remote_local (bool):   keep local copies of remote files (default False)
        nodeps (bool):              ignore dependencies (default False)
        keep_target_files (bool):   do not adjust the paths of given target files relative to the working directory.
        allowed_rules (set):        restrict allowed rules to the given set. If None or empty, all rules are used.
        jobscript (str):            path to a custom shell script template for cluster jobs (default None)
        greediness (float):         set the greediness of scheduling. This value between 0 and 1 determines how careful jobs are selected for execution. The default value (0.5 if prioritytargets are used, 1.0 else) provides the best speed and still acceptable scheduling quality.
        overwrite_shellcmd (str):   a shell command that shall be executed instead of those given in the workflow. This is for debugging purposes only.
        updated_files(list):        a list that will be filled with the files that are updated or created during the workflow execution
        verbose (bool):             show additional debug output (default False)
        max_jobs_per_second (int):  maximal number of cluster/drmaa jobs per second, None to impose no limit (default None)
        restart_times (int):        number of times to restart failing jobs (default 0)
        attempt (int):              initial value of Job.attempt. This is intended for internal use only (default 1).
        force_use_threads:          whether to force use of threads over processes. helpful if shared memory is full or unavailable (default False)
        use_conda (bool):           use conda environments for each job (defined with conda directive of rules)
        use_singularity (bool):     run jobs in singularity containers (if defined with singularity directive)
        use_env_modules (bool):     load environment modules if defined in rules
        singularity_args (str):     additional arguments to pass to singularity
        conda_prefix (str):         the directory in which conda environments will be created (default None)
        conda_cleanup_pkgs (snakemake.deployment.conda.CondaCleanupMode):
                                    whether to clean up conda tarballs after env creation (default None), valid values: "tarballs", "cache"
        singularity_prefix (str):   the directory to which singularity images will be pulled (default None)
        shadow_prefix (str):        prefix for shadow directories. The job-specific shadow directories will be created in $SHADOW_PREFIX/shadow/ (default None)
        wms-monitor (str):          workflow management system monitor. Send post requests to the specified (server/IP). (default None)
        conda_create_envs_only (bool):    if specified, only builds the conda environments specified for each job, then exits.
        list_conda_envs (bool):     list conda environments and their location on disk.
        mode (snakemake.common.Mode): execution mode
        wrapper_prefix (str):       prefix for wrapper script URLs (default None)
        kubernetes (str):           submit jobs to kubernetes, using the given namespace.
        container_image (str):      Docker image to use, e.g., for kubernetes.
        default_remote_provider (str): default remote provider to use instead of local files (e.g. S3, GS)
        default_remote_prefix (str): prefix for default remote provider (e.g. name of the bucket).
        tibanna (bool):             submit jobs to AWS cloud using Tibanna.
        tibanna_sfn (str):          Step function (Unicorn) name of Tibanna (e.g. tibanna_unicorn_monty). This must be deployed first using tibanna cli.
        google_lifesciences (bool): submit jobs to Google Cloud Life Sciences (pipelines API).
        google_lifesciences_regions (list): a list of regions (e.g., us-east1)
        google_lifesciences_location (str): Life Sciences API location (e.g., us-central1)
        google_lifesciences_cache (bool): save a cache of the compressed working directories in Google Cloud Storage for later usage.
        tes (str):                  Execute workflow tasks on GA4GH TES server given by url.
        precommand (str):           commands to run on AWS cloud before the snakemake command (e.g. wget, git clone, unzip, etc). Use with --tibanna.
        preemption_default (int):   set a default number of preemptible instance retries (for Google Life Sciences executor only)
        preemptible_rules (list):    define custom preemptible instance retries for specific rules (for Google Life Sciences executor only)
        tibanna_config (list):      Additional tibanna config e.g. --tibanna-config spot_instance=true subnet=<subnet_id> security group=<security_group_id>
        assume_shared_fs (bool):    assume that cluster nodes share a common filesystem (default true).
        cluster_status (str):       status command for cluster execution. If None, Snakemake will rely on flag files. Otherwise, it expects the command to return "success", "failure" or "running" when executing with a cluster jobid as single argument.
        export_cwl (str):           Compile workflow to CWL and save to given file
        log_handler (function):     redirect snakemake output to this custom log handler, a function that takes a log message dictionary (see below) as its only argument (default None). The log message dictionary for the log handler has to following entries:
        keep_incomplete (bool):     keep incomplete output files of failed jobs
        edit_notebook (object):     "notebook.Listen" object to configuring notebook server for interactive editing of a rule notebook. If None, do not edit.
        scheduler (str):            Select scheduling algorithm (default ilp)
        scheduler_ilp_solver (str): Set solver for ilp scheduler.
        overwrite_groups (dict):    Rule to group assignments (default None)
        group_components (dict):    Number of connected components given groups shall span before being split up (1 by default if empty)
        log_handler (list):         redirect snakemake output to this list of custom log handler, each a function that takes a log message dictionary (see below) as its only argument (default []). The log message dictionary for the log handler has to following entries:

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
    assert not immediate_submit or (
        immediate_submit and notemp
    ), "immediate_submit has to be combined with notemp (it does not support temp file handling)"

    if tibanna:
        assume_shared_fs = False
        default_remote_provider = "S3"
        default_remote_prefix = default_remote_prefix.rstrip("/")
        assert (
            default_remote_prefix
        ), "default_remote_prefix needed if tibanna is specified"
        assert tibanna_sfn, "tibanna_sfn needed if tibanna is specified"
        if tibanna_config:
            tibanna_config_dict = dict()
            for cf in tibanna_config:
                k, v = cf.split("=")
                if v == "true":
                    v = True
                elif v == "false":
                    v = False
                elif v.isnumeric():
                    v = int(v)
                else:
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                tibanna_config_dict.update({k: v})
            tibanna_config = tibanna_config_dict

    # Google Cloud Life Sciences API uses compute engine and storage
    if google_lifesciences:
        assume_shared_fs = False
        default_remote_provider = "GS"
        default_remote_prefix = default_remote_prefix.rstrip("/")

    # Currently preemptible instances only supported for Google LifeSciences Executor
    if preemption_default or preemptible_rules and not google_lifesciences:
        logger.warning(
            "Preemptible instances are only available for the Google Life Sciences Executor."
        )

    if updated_files is None:
        updated_files = list()

    if (
        cluster
        or cluster_sync
        or drmaa
        or tibanna
        or kubernetes
        or google_lifesciences
        or tes
    ):
        cores = None
    else:
        nodes = None

    if isinstance(cluster_config, str):
        # Loading configuration from one file is still supported for
        # backward compatibility
        cluster_config = [cluster_config]
    if cluster_config:
        # Load all configuration files
        configs = [load_configfile(f) for f in cluster_config]
        # Merge in the order as specified, overriding earlier values with
        # later ones
        cluster_config_content = configs[0]
        for other in configs[1:]:
            update_config(cluster_config_content, other)
    else:
        cluster_config_content = dict()

    run_local = not (
        cluster
        or cluster_sync
        or drmaa
        or kubernetes
        or tibanna
        or google_lifesciences
        or tes
    )
    if run_local:
        if not dryrun:
            # clean up all previously recorded jobids.
            shell.cleanup()
    else:
        if edit_notebook:
            raise WorkflowError(
                "Notebook edit mode is only allowed with local execution."
            )

    # force thread use for any kind of cluster
    use_threads = (
        force_use_threads
        or (os.name not in ["posix", "nt"])
        or cluster
        or cluster_sync
        or drmaa
    )

    if not keep_logger:
        stdout = (
            (
                dryrun
                and not (printdag or printd3dag or printrulegraph or printfilegraph)
            )
            or listrules
            or list_target_rules
            or list_resources
        )

        setup_logger(
            handler=log_handler,
            quiet=quiet,
            printreason=printreason,
            printshellcmds=printshellcmds,
            debug_dag=debug_dag,
            nocolor=nocolor,
            stdout=stdout,
            debug=verbose,
            use_threads=use_threads,
            mode=mode,
            show_failed_logs=show_failed_logs,
            wms_monitor=wms_monitor,
        )

    if greediness is None:
        greediness = 0.5 if prioritytargets else 1.0
    else:
        if not (0 <= greediness <= 1.0):
            logger.error("Error: greediness must be a float between 0 and 1.")
            return False

    if not os.path.exists(snakefile):
        logger.error('Error: Snakefile "{}" not found.'.format(snakefile))
        return False
    snakefile = os.path.abspath(snakefile)

    cluster_mode = (
        (cluster is not None) + (cluster_sync is not None) + (drmaa is not None)
    )
    if cluster_mode > 1:
        logger.error("Error: cluster and drmaa args are mutually exclusive")
        return False

    if debug and (cluster_mode or cores is not None and cores > 1):
        logger.error(
            "Error: debug mode cannot be used with more than one core or cluster execution."
        )
        return False

    overwrite_config = dict()
    if configfiles is None:
        configfiles = []
    for f in configfiles:
        # get values to override. Later configfiles override earlier ones.
        overwrite_config.update(load_configfile(f))
    # convert provided paths to absolute paths
    configfiles = list(map(os.path.abspath, configfiles))

    # directly specified elements override any configfiles
    if config:
        overwrite_config.update(config)
        if config_args is None:
            config_args = unparse_config(config)

    if workdir:
        olddir = os.getcwd()
        if not os.path.exists(workdir):
            logger.info("Creating specified working directory {}.".format(workdir))
            os.makedirs(workdir)
        workdir = os.path.abspath(workdir)
        os.chdir(workdir)

    logger.setup_logfile()

    try:
        # handle default remote provider
        _default_remote_provider = None
        if default_remote_provider is not None:
            try:
                rmt = importlib.import_module(
                    "snakemake.remote." + default_remote_provider
                )
            except ImportError as e:
                raise WorkflowError("Unknown default remote provider.")
            if rmt.RemoteProvider.supports_default:
                _default_remote_provider = rmt.RemoteProvider(
                    keep_local=True, is_default=True
                )
            else:
                raise WorkflowError(
                    "Remote provider {} does not (yet) support to "
                    "be used as default provider."
                )

        workflow = Workflow(
            snakefile=snakefile,
            jobscript=jobscript,
            overwrite_shellcmd=overwrite_shellcmd,
            overwrite_config=overwrite_config,
            overwrite_workdir=workdir,
            overwrite_configfiles=configfiles,
            overwrite_clusterconfig=cluster_config_content,
            overwrite_threads=overwrite_threads,
            overwrite_scatter=overwrite_scatter,
            overwrite_groups=overwrite_groups,
            group_components=group_components,
            config_args=config_args,
            debug=debug,
            verbose=verbose,
            use_conda=use_conda or list_conda_envs or conda_cleanup_envs,
            use_singularity=use_singularity,
            use_env_modules=use_env_modules,
            conda_frontend=conda_frontend,
            conda_prefix=conda_prefix,
            conda_cleanup_pkgs=conda_cleanup_pkgs,
            singularity_prefix=singularity_prefix,
            shadow_prefix=shadow_prefix,
            singularity_args=singularity_args,
            scheduler_type=scheduler,
            scheduler_ilp_solver=scheduler_ilp_solver,
            mode=mode,
            wrapper_prefix=wrapper_prefix,
            printshellcmds=printshellcmds,
            restart_times=restart_times,
            attempt=attempt,
            default_remote_provider=_default_remote_provider,
            default_remote_prefix=default_remote_prefix,
            run_local=run_local,
            default_resources=default_resources,
            cache=cache,
            cores=cores,
            nodes=nodes,
            resources=resources,
            edit_notebook=edit_notebook,
            envvars=envvars,
            max_inventory_wait_time=max_inventory_wait_time,
        )
        success = True

        workflow.include(
            snakefile, overwrite_first_rule=True, print_compilation=print_compilation
        )
        workflow.check()

        if not print_compilation:
            if lint:
                success = not workflow.lint(json=lint == "json")
            elif listrules:
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
                    local_cores=local_cores,
                    cache=cache,
                    overwrite_threads=overwrite_threads,
                    overwrite_scatter=overwrite_scatter,
                    default_resources=default_resources,
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
                    restart_times=restart_times,
                    attempt=attempt,
                    lock=lock,
                    unlock=unlock,
                    cleanup_metadata=cleanup_metadata,
                    conda_cleanup_envs=conda_cleanup_envs,
                    cleanup_shadow=cleanup_shadow,
                    cleanup_scripts=cleanup_scripts,
                    force_incomplete=force_incomplete,
                    ignore_incomplete=ignore_incomplete,
                    latency_wait=latency_wait,
                    verbose=verbose,
                    notemp=notemp,
                    keep_remote_local=keep_remote_local,
                    nodeps=nodeps,
                    jobscript=jobscript,
                    greediness=greediness,
                    no_hooks=no_hooks,
                    overwrite_shellcmd=overwrite_shellcmd,
                    config=config,
                    config_args=config_args,
                    cluster_config=cluster_config,
                    keep_logger=True,
                    force_use_threads=use_threads,
                    use_conda=use_conda,
                    use_singularity=use_singularity,
                    use_env_modules=use_env_modules,
                    conda_prefix=conda_prefix,
                    conda_cleanup_pkgs=conda_cleanup_pkgs,
                    singularity_prefix=singularity_prefix,
                    shadow_prefix=shadow_prefix,
                    singularity_args=singularity_args,
                    scheduler=scheduler,
                    scheduler_ilp_solver=scheduler_ilp_solver,
                    list_conda_envs=list_conda_envs,
                    kubernetes=kubernetes,
                    container_image=container_image,
                    conda_create_envs_only=conda_create_envs_only,
                    default_remote_provider=default_remote_provider,
                    default_remote_prefix=default_remote_prefix,
                    tibanna=tibanna,
                    tibanna_sfn=tibanna_sfn,
                    google_lifesciences=google_lifesciences,
                    google_lifesciences_regions=google_lifesciences_regions,
                    google_lifesciences_location=google_lifesciences_location,
                    google_lifesciences_cache=google_lifesciences_cache,
                    tes=tes,
                    precommand=precommand,
                    preemption_default=preemption_default,
                    preemptible_rules=preemptible_rules,
                    tibanna_config=tibanna_config,
                    assume_shared_fs=assume_shared_fs,
                    cluster_status=cluster_status,
                    max_jobs_per_second=max_jobs_per_second,
                    max_status_checks_per_second=max_status_checks_per_second,
                    overwrite_groups=overwrite_groups,
                    group_components=group_components,
                    max_inventory_wait_time=max_inventory_wait_time,
                )
                success = workflow.execute(
                    targets=targets,
                    dryrun=dryrun,
                    generate_unit_tests=generate_unit_tests,
                    touch=touch,
                    scheduler_type=scheduler,
                    scheduler_ilp_solver=scheduler_ilp_solver,
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
                    printfilegraph=printfilegraph,
                    printdag=printdag,
                    cluster=cluster,
                    cluster_sync=cluster_sync,
                    jobname=jobname,
                    drmaa=drmaa,
                    drmaa_log_dir=drmaa_log_dir,
                    kubernetes=kubernetes,
                    container_image=container_image,
                    tibanna=tibanna,
                    tibanna_sfn=tibanna_sfn,
                    google_lifesciences=google_lifesciences,
                    google_lifesciences_regions=google_lifesciences_regions,
                    google_lifesciences_location=google_lifesciences_location,
                    google_lifesciences_cache=google_lifesciences_cache,
                    tes=tes,
                    precommand=precommand,
                    preemption_default=preemption_default,
                    preemptible_rules=preemptible_rules,
                    tibanna_config=tibanna_config,
                    max_jobs_per_second=max_jobs_per_second,
                    max_status_checks_per_second=max_status_checks_per_second,
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
                    list_untracked=list_untracked,
                    list_conda_envs=list_conda_envs,
                    summary=summary,
                    archive=archive,
                    delete_all_output=delete_all_output,
                    delete_temp_output=delete_temp_output,
                    latency_wait=latency_wait,
                    wait_for_files=wait_for_files,
                    detailed_summary=detailed_summary,
                    nolock=not lock,
                    unlock=unlock,
                    notemp=notemp,
                    keep_remote_local=keep_remote_local,
                    nodeps=nodeps,
                    keep_target_files=keep_target_files,
                    cleanup_metadata=cleanup_metadata,
                    conda_cleanup_envs=conda_cleanup_envs,
                    cleanup_shadow=cleanup_shadow,
                    cleanup_scripts=cleanup_scripts,
                    subsnakemake=subsnakemake,
                    updated_files=updated_files,
                    allowed_rules=allowed_rules,
                    greediness=greediness,
                    no_hooks=no_hooks,
                    force_use_threads=use_threads,
                    conda_create_envs_only=conda_create_envs_only,
                    assume_shared_fs=assume_shared_fs,
                    cluster_status=cluster_status,
                    report=report,
                    report_stylesheet=report_stylesheet,
                    export_cwl=export_cwl,
                    batch=batch,
                    keepincomplete=keep_incomplete,
                    keepmetadata=keep_metadata,
                    executesubworkflows=execute_subworkflows,
                )

    except BrokenPipeError:
        # ignore this exception and stop. It occurs if snakemake output is piped into less and less quits before reading the whole output.
        # in such a case, snakemake shall stop scheduling and quit with error 1
        success = False
    except (Exception, BaseException) as ex:
        if "workflow" in locals():
            print_exception(ex, workflow.linemaps)
        else:
            print_exception(ex, dict())
        success = False

    if workdir:
        os.chdir(olddir)
    if "workflow" in locals() and workflow.persistence:
        workflow.persistence.unlock()
    if not keep_logger:
        logger.cleanup()
    return success


def parse_set_threads(args):
    return parse_set_ints(
        args.set_threads,
        "Invalid threads definition: entries have to be defined as RULE=THREADS pairs "
        "(with THREADS being a positive integer).",
    )


def parse_set_scatter(args):
    return parse_set_ints(
        args.set_scatter,
        "Invalid scatter definition: entries have to be defined as NAME=SCATTERITEMS pairs "
        "(with SCATTERITEMS being a positive integer).",
    )


def parse_set_ints(arg, errmsg):
    assignments = dict()
    if arg is not None:
        for entry in arg:
            key, value = parse_key_value_arg(entry, errmsg=errmsg)
            try:
                value = int(value)
            except ValueError:
                raise ValueError(errmsg)
            if value < 0:
                raise ValueError(errmsg)
            assignments[key] = value
    return assignments


def parse_batch(args):
    errmsg = "Invalid batch definition: batch entry has to be defined as RULE=BATCH/BATCHES (with integers BATCH <= BATCHES, BATCH >= 1)."
    if args.batch is not None:
        rule, batchdef = parse_key_value_arg(args.batch, errmsg=errmsg)
        try:
            batch, batches = batchdef.split("/")
            batch = int(batch)
            batches = int(batches)
        except ValueError:
            raise ValueError(errmsg)
        if batch > batches or batch < 1:
            raise ValueError(errmsg)
        return Batch(rule, batch, batches)
    return None


def parse_groups(args):
    errmsg = "Invalid groups definition: entries have to be defined as RULE=GROUP pairs"
    overwrite_groups = dict()
    if args.groups is not None:
        for entry in args.groups:
            rule, group = parse_key_value_arg(entry, errmsg=errmsg)
            overwrite_groups[rule] = group
    return overwrite_groups


def parse_group_components(args):
    errmsg = "Invalid group components definition: entries have to be defined as GROUP=COMPONENTS pairs (with COMPONENTS being a positive integer)"
    group_components = dict()
    if args.group_components is not None:
        for entry in args.group_components:
            group, count = parse_key_value_arg(entry, errmsg=errmsg)
            try:
                count = int(count)
            except ValueError:
                raise ValueError(errmsg)
            if count <= 0:
                raise ValueError(errmsg)
            group_components[group] = count
    return group_components


def parse_key_value_arg(arg, errmsg):
    try:
        key, val = arg.split("=", 1)
    except ValueError:
        raise ValueError(errmsg + " Unparseable value: %r." % arg)
    return key, val


def _bool_parser(value):
    if value == "True":
        return True
    elif value == "False":
        return False
    raise ValueError


def parse_config(args):
    """Parse config from args."""
    import yaml

    yaml_base_load = lambda s: yaml.load(s, Loader=yaml.loader.BaseLoader)
    parsers = [int, float, _bool_parser, yaml_base_load, str]
    config = dict()
    if args.config is not None:
        valid = re.compile(r"[a-zA-Z_]\w*$")
        for entry in args.config:
            key, val = parse_key_value_arg(
                entry,
                errmsg="Invalid config definition: Config entries have to be defined as name=value pairs.",
            )
            if not valid.match(key):
                raise ValueError(
                    "Invalid config definition: Config entry must start with a valid identifier."
                )
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


def unparse_config(config):
    if not isinstance(config, dict):
        raise ValueError("config is not a dict")
    items = []
    for key, value in config.items():
        if isinstance(value, dict):
            raise ValueError("config may only be a flat dict")
        encoded = "'{}'".format(value) if isinstance(value, str) else value
        items.append("{}={}".format(key, encoded))
    return items


APPDIRS = None


def get_appdirs():
    global APPDIRS
    if APPDIRS is None:
        from appdirs import AppDirs

        APPDIRS = AppDirs("snakemake", "snakemake")
    return APPDIRS


def get_profile_file(profile, file, return_default=False):
    dirs = get_appdirs()
    if os.path.isabs(profile):
        search_dirs = [os.path.dirname(profile)]
        profile = os.path.basename(profile)
    else:
        search_dirs = [os.getcwd(), dirs.user_config_dir, dirs.site_config_dir]
    get_path = lambda d: os.path.join(d, profile, file)
    for d in search_dirs:
        p = get_path(d)
        if os.path.exists(p):
            return p

    if return_default:
        return file
    return None


def get_argument_parser(profile=None):
    """Generate and return argument parser."""
    import configargparse
    from configargparse import YAMLConfigFileParser

    dirs = get_appdirs()
    config_files = []
    if profile:
        if profile == "":
            print("Error: invalid profile name.", file=sys.stderr)
            exit(1)

        config_file = get_profile_file(profile, "config.yaml")
        if config_file is None:
            print(
                "Error: profile given but no config.yaml found. "
                "Profile has to be given as either absolute path, relative "
                "path or name of a directory available in either "
                "{site} or {user}.".format(
                    site=dirs.site_config_dir, user=dirs.user_config_dir
                ),
                file=sys.stderr,
            )
            exit(1)
        config_files = [config_file]

    parser = configargparse.ArgumentParser(
        description="Snakemake is a Python based language and execution "
        "environment for GNU Make-like workflows.",
        formatter_class=ArgumentDefaultsHelpFormatter,
        default_config_files=config_files,
        config_file_parser_class=YAMLConfigFileParser,
    )

    group_exec = parser.add_argument_group("EXECUTION")

    group_exec.add_argument(
        "target",
        nargs="*",
        default=None,
        help="Targets to build. May be rules or files.",
    )

    group_exec.add_argument(
        "--dry-run",
        "--dryrun",
        "-n",
        dest="dryrun",
        action="store_true",
        help="Do not execute anything, and display what would be done. "
        "If you have a very large workflow, use --dry-run --quiet to just "
        "print a summary of the DAG of jobs.",
    )

    group_exec.add_argument(
        "--profile",
        help="""
                        Name of profile to use for configuring
                        Snakemake. Snakemake will search for a corresponding
                        folder in {} and {}. Alternatively, this can be an
                        absolute or relative path.
                        The profile folder has to contain a file 'config.yaml'.
                        This file can be used to set default values for command
                        line options in YAML format. For example,
                        '--cluster qsub' becomes 'cluster: qsub' in the YAML
                        file. Profiles can be obtained from
                        https://github.com/snakemake-profiles.
                        """.format(
            dirs.site_config_dir, dirs.user_config_dir
        ),
    )

    group_exec.add_argument(
        "--cache",
        nargs="*",
        metavar="RULE",
        help="Store output files of given rules in a central cache given by the environment "
        "variable $SNAKEMAKE_OUTPUT_CACHE. Likewise, retrieve output files of the given rules "
        "from this cache if they have been created before (by anybody writing to the same cache), "
        "instead of actually executing the rules. Output files are identified by hashing all "
        "steps, parameters and software stack (conda envs or containers) needed to create them.",
    )

    group_exec.add_argument(
        "--snakefile",
        "-s",
        metavar="FILE",
        help=(
            "The workflow definition in form of a snakefile."
            "Usually, you should not need to specify this. "
            "By default, Snakemake will search for {} "
            "beneath the current working "
            "directory, in this order. "
            "Only if you definitely want a different layout, "
            "you need to use this parameter."
        ).format(", ".join(map("'{}'".format, SNAKEFILE_CHOICES))),
    )
    group_exec.add_argument(
        "--cores",
        "--jobs",
        "-j",
        action="store",
        const=available_cpu_count(),
        nargs="?",
        metavar="N",
        help=(
            "Use at most N CPU cores/jobs in parallel. "
            "If N is omitted or 'all', the limit is set to the number of "
            "available CPU cores."
        ),
    )
    group_exec.add_argument(
        "--local-cores",
        action="store",
        default=available_cpu_count(),
        metavar="N",
        type=int,
        help=(
            "In cluster mode, use at most N cores of the host machine in parallel "
            " (default: number of CPU cores of the host). The cores are used to execute "
            "local rules. This option is ignored when not in cluster mode."
        ),
    )
    group_exec.add_argument(
        "--resources",
        "--res",
        nargs="*",
        metavar="NAME=INT",
        help=(
            "Define additional resources that shall constrain the scheduling "
            "analogously to threads (see above). A resource is defined as "
            "a name and an integer value. E.g. --resources mem_mb=1000. Rules can "
            "use resources by defining the resource keyword, e.g. "
            "resources: mem_mb=600. If now two rules require 600 of the resource "
            "'mem_mb' they won't be run in parallel by the scheduler."
        ),
    )
    group_exec.add_argument(
        "--set-threads",
        metavar="RULE=THREADS",
        nargs="+",
        help="Overwrite thread usage of rules. This allows to fine-tune workflow "
        "parallelization. In particular, this is helpful to target certain cluster nodes "
        "by e.g. shifting a rule to use more, or less threads than defined in the workflow. "
        "Thereby, THREADS has to be a positive integer, and RULE has to be the name of the rule.",
    )
    group_exec.add_argument(
        "--set-scatter",
        metavar="NAME=SCATTERITEMS",
        nargs="+",
        help="Overwrite number of scatter items of scattergather processes. This allows to fine-tune "
        "workflow parallelization. Thereby, SCATTERITEMS has to be a positive integer, and NAME has to be "
        "the name of the scattergather process defined via a scattergather directive in the workflow.",
    )
    group_exec.add_argument(
        "--default-resources",
        "--default-res",
        nargs="*",
        metavar="NAME=INT",
        help=(
            "Define default values of resources for rules that do not define their own values. "
            "In addition to plain integers, python expressions over inputsize are allowed (e.g. '2*input.size_mb')."
            "When specifying this without any arguments (--default-resources), it defines 'mem_mb=max(2*input.size_mb, 1000)' "
            "'disk_mb=max(2*input.size_mb, 1000)', i.e., default disk and mem usage is twice the input file size but at least 1GB."
        ),
    )

    group_exec.add_argument(
        "--preemption-default",
        type=int,
        default=None,
        help=(
            "A preemptible instance can be requested when using the Google Life Sciences API. If you set a --preemption-default,"
            "all rules will be subject to the default. Specifically, this integer is the number of restart attempts that will be "
            "made given that the instance is killed unexpectedly. Note that preemptible instances have a maximum running time of 24 "
            "hours. If you want to set preemptible instances for only a subset of rules, use --preemptible-rules instead."
        ),
    )

    group_exec.add_argument(
        "--preemptible-rules",
        nargs="+",
        default=None,
        help=(
            "A preemptible instance can be requested when using the Google Life Sciences API. If you want to use these instances "
            "for a subset of your rules, you can use --preemptible-rules and then specify a list of rule and integer pairs, where "
            "each integer indicates the number of restarts to use for the rule's instance in the case that the instance is "
            "terminated unexpectedly. --preemptible-rules can be used in combination with --preemption-default, and will take "
            "priority. Note that preemptible instances have a maximum running time of 24. If you want to apply a consistent "
            "number of retries across all your rules, use --premption-default instead. "
            "Example: snakemake --preemption-default 10 --preemptible-rules map_reads=3 call_variants=0"
        ),
    )

    group_exec.add_argument(
        "--config",
        "-C",
        nargs="*",
        metavar="KEY=VALUE",
        help=(
            "Set or overwrite values in the workflow config object. "
            "The workflow config object is accessible as variable config inside "
            "the workflow. Default values can be set by providing a JSON file "
            "(see Documentation)."
        ),
    )
    group_exec.add_argument(
        "--configfile",
        "--configfiles",
        nargs="+",
        metavar="FILE",
        help=(
            "Specify or overwrite the config file of the workflow (see the docs). "
            "Values specified in JSON or YAML format are available in the global config "
            "dictionary inside the workflow. Multiple files overwrite each other in "
            "the given order."
        ),
    )
    group_exec.add_argument(
        "--envvars",
        nargs="+",
        metavar="VARNAME",
        help="Environment variables to pass to cloud jobs.",
    )
    group_exec.add_argument(
        "--directory",
        "-d",
        metavar="DIR",
        action="store",
        help=(
            "Specify working directory (relative paths in "
            "the snakefile will use this as their origin)."
        ),
    )
    group_exec.add_argument(
        "--touch",
        "-t",
        action="store_true",
        help=(
            "Touch output files (mark them up to date without really "
            "changing them) instead of running their commands. This is "
            "used to pretend that the rules were executed, in order to "
            "fool future invocations of snakemake. Fails if a file does "
            "not yet exist. Note that this will only touch files that would "
            "otherwise be recreated by Snakemake (e.g. because their input "
            "files are newer). For enforcing a touch, combine this with "
            "--force, --forceall, or --forcerun. Note however that you loose "
            "the provenance information when the files have been created in "
            "realitiy. Hence, this should be used only as a last resort."
        ),
    )
    group_exec.add_argument(
        "--keep-going",
        "-k",
        action="store_true",
        help="Go on with independent jobs if a job fails.",
    )
    group_exec.add_argument(
        "--force",
        "-f",
        action="store_true",
        help=(
            "Force the execution of the selected target or the first rule "
            "regardless of already created output."
        ),
    )
    group_exec.add_argument(
        "--forceall",
        "-F",
        action="store_true",
        help=(
            "Force the execution of the selected (or the first) rule and "
            "all rules it is dependent on regardless of already created "
            "output."
        ),
    )
    group_exec.add_argument(
        "--forcerun",
        "-R",
        nargs="*",
        metavar="TARGET",
        help=(
            "Force the re-execution or creation of the given rules or files."
            " Use this option if you changed a rule and want to have all its "
            "output in your workflow updated."
        ),
    )
    group_exec.add_argument(
        "--prioritize",
        "-P",
        nargs="+",
        metavar="TARGET",
        help=(
            "Tell the scheduler to assign creation of given targets "
            "(and all their dependencies) highest priority. (EXPERIMENTAL)"
        ),
    )
    group_exec.add_argument(
        "--batch",
        metavar="RULE=BATCH/BATCHES",
        help=(
            "Only create the given BATCH of the input files of the given RULE. "
            "This can be used to iteratively run parts of very large workflows. "
            "Only the execution plan of the relevant part of the workflow has to "
            "be calculated, thereby speeding up DAG computation. "
            "It is recommended to provide the most suitable rule for batching when "
            "documenting a workflow. It should be some aggregating rule that "
            "would be executed only once, and has a large number of input files. "
            "For example, it can be a rule that aggregates over samples."
        ),
    )
    group_exec.add_argument(
        "--until",
        "-U",
        nargs="+",
        metavar="TARGET",
        help=(
            "Runs the pipeline until it reaches the specified rules or "
            "files. Only runs jobs that are dependencies of the specified "
            "rule or files, does not run sibling DAGs. "
        ),
    )
    group_exec.add_argument(
        "--omit-from",
        "-O",
        nargs="+",
        metavar="TARGET",
        help=(
            "Prevent the execution or creation of the given rules or files "
            "as well as any rules or files that are downstream of these targets "
            "in the DAG. Also runs jobs in sibling DAGs that are independent of the "
            "rules or files specified here."
        ),
    )
    group_exec.add_argument(
        "--rerun-incomplete",
        "--ri",
        action="store_true",
        help=("Re-run all " "jobs the output of which is recognized as incomplete."),
    )
    group_exec.add_argument(
        "--shadow-prefix",
        metavar="DIR",
        help=(
            "Specify a directory in which the 'shadow' directory is created. "
            "If not supplied, the value is set to the '.snakemake' directory relative "
            "to the working directory."
        ),
    )

    try:
        import pulp

        lp_solvers = pulp.list_solvers(onlyAvailable=True)
    except ImportError:
        # Dummy list for the case that pulp is not available
        # This only happend when building docs.
        lp_solvers = ["COIN_CMD"]
    recommended_lp_solver = "COIN_CMD"

    group_exec.add_argument(
        "--scheduler",
        default="greedy" if recommended_lp_solver not in lp_solvers else "ilp",
        nargs="?",
        choices=["ilp", "greedy"],
        help=(
            "Specifies if jobs are selected by a greedy algorithm or by solving an ilp. "
            "The ilp scheduler aims to reduce runtime and hdd usage by best possible use of resources."
        ),
    )
    group_exec.add_argument(
        "--wms-monitor",
        action="store",
        nargs="?",
        help=(
            "IP and port of workflow management system to monitor the execution of snakemake (e.g. http://127.0.0.1:5000"
        ),
    )

    group_exec.add_argument(
        "--scheduler-ilp-solver",
        default=recommended_lp_solver,
        choices=lp_solvers,
        help=("Specifies solver to be utilized when selecting ilp-scheduler."),
    )

    group_exec.add_argument(
        "--no-subworkflows",
        "--nosw",
        action="store_true",
        help=("Do not evaluate or execute subworkflows."),
    )

    # TODO add group_partitioning, allowing to define --group rulename=groupname.
    # i.e. setting groups via the CLI for improving cluster performance given
    # available resources.
    # TODO add an additional flag --group-components groupname=3, allowing to set the
    # number of connected components a group is allowed to span. By default, this is 1
    # (as now), but the flag allows to extend this. This can be used to run e.g.
    # 3 jobs of the same rule in the same group, although they are not connected.
    # Can be helpful for putting together many small jobs or benefitting of shared memory
    # setups.

    group_group = parser.add_argument_group("GROUPING")
    group_group.add_argument(
        "--groups",
        nargs="+",
        help="Assign rules to groups (this overwrites any "
        "group definitions from the workflow).",
    )
    group_group.add_argument(
        "--group-components",
        nargs="+",
        help="Set the number of connected components a group is "
        "allowed to span. By default, this is 1, but this flag "
        "allows to extend this. This can be used to run e.g. 3 "
        "jobs of the same rule in the same group, although they "
        "are not connected. It can be helpful for putting together "
        "many small jobs or benefitting of shared memory setups.",
    )

    group_report = parser.add_argument_group("REPORTS")

    group_report.add_argument(
        "--report",
        nargs="?",
        const="report.html",
        metavar="FILE",
        help="Create an HTML report with results and statistics. "
        "This can be either a .html file or a .zip file. "
        "In the former case, all results are embedded into the .html (this only works for small data). "
        "In the latter case, results are stored along with a file report.html in the zip archive. "
        "If no filename is given, an embedded report.html is the default.",
    )
    group_report.add_argument(
        "--report-stylesheet",
        metavar="CSSFILE",
        help="Custom stylesheet to use for report. In particular, this can be used for "
        "branding the report with e.g. a custom logo, see docs.",
    )

    group_notebooks = parser.add_argument_group("NOTEBOOKS")

    group_notebooks.add_argument(
        "--edit-notebook",
        metavar="TARGET",
        help="Interactively edit the notebook associated with the rule used to generate the given target file. "
        "This will start a local jupyter notebook server. "
        "Any changes to the notebook should be saved, and the server has to be stopped by "
        "closing the notebook and hitting the 'Quit' button on the jupyter dashboard. "
        "Afterwards, the updated notebook will be automatically stored in the path defined in the rule. "
        "If the notebook is not yet present, this will create an empty draft. ",
    )
    group_notebooks.add_argument(
        "--notebook-listen",
        metavar="IP:PORT",
        default="localhost:8888",
        help="The IP address and PORT the notebook server used for editing the notebook (--edit-notebook) will listen on.",
    )

    group_utils = parser.add_argument_group("UTILITIES")
    group_utils.add_argument(
        "--lint",
        nargs="?",
        const="text",
        choices=["text", "json"],
        help="Perform linting on the given workflow. This will print snakemake "
        "specific suggestions to improve code quality (work in progress, more lints "
        "to be added in the future). If no argument is provided, plain text output is used.",
    )
    group_utils.add_argument(
        "--generate-unit-tests",
        nargs="?",
        const=".tests/unit",
        metavar="TESTPATH",
        help="Automatically generate unit tests for each workflow rule. "
        "This assumes that all input files of each job are already present. "
        "Rules without a job with present input files will be skipped (a warning will be issued). "
        "For each rule, one test case will be "
        "created in the specified test folder (.tests/unit by default). After "
        "successfull execution, tests can be run with "
        "'pytest TESTPATH'.",
    )
    group_utils.add_argument(
        "--export-cwl",
        action="store",
        metavar="FILE",
        help="Compile workflow to CWL and store it in given FILE.",
    )
    group_utils.add_argument(
        "--list",
        "-l",
        action="store_true",
        help="Show available rules in given Snakefile.",
    )
    group_utils.add_argument(
        "--list-target-rules",
        "--lt",
        action="store_true",
        help="Show available target rules in given Snakefile.",
    )
    group_utils.add_argument(
        "--dag",
        action="store_true",
        help="Do not execute anything and print the directed "
        "acyclic graph of jobs in the dot language. Recommended "
        "use on Unix systems: snakemake --dag | dot | display"
        "Note print statements in your Snakefile may interfere "
        "with visualization.",
    )
    group_utils.add_argument(
        "--rulegraph",
        action="store_true",
        help="Do not execute anything and print the dependency graph "
        "of rules in the dot language. This will be less "
        "crowded than above DAG of jobs, but also show less information. "
        "Note that each rule is displayed once, hence the displayed graph will be "
        "cyclic if a rule appears in several steps of the workflow. "
        "Use this if above option leads to a DAG that is too large. "
        "Recommended use on Unix systems: snakemake --rulegraph | dot | display"
        "Note print statements in your Snakefile may interfere "
        "with visualization.",
    )
    group_utils.add_argument(
        "--filegraph",
        action="store_true",
        help="Do not execute anything and print the dependency graph "
        "of rules with their input and output files in the dot language. "
        "This is an intermediate solution between above DAG of jobs and the rule graph. "
        "Note that each rule is displayed once, hence the displayed graph will be "
        "cyclic if a rule appears in several steps of the workflow. "
        "Use this if above option leads to a DAG that is too large. "
        "Recommended use on Unix systems: snakemake --filegraph | dot | display"
        "Note print statements in your Snakefile may interfere "
        "with visualization.",
    )
    group_utils.add_argument(
        "--d3dag",
        action="store_true",
        help="Print the DAG in D3.js compatible JSON format.",
    )
    group_utils.add_argument(
        "--summary",
        "-S",
        action="store_true",
        help="Print a summary of all files created by the workflow. The "
        "has the following columns: filename, modification time, "
        "rule version, status, plan.\n"
        "Thereby rule version contains the version"
        "the file was created with (see the version keyword of rules), and "
        "status denotes whether the file is missing, its input files are "
        "newer or if version or implementation of the rule changed since "
        "file creation. Finally the last column denotes whether the file "
        "will be updated or created during the next workflow execution.",
    )
    group_utils.add_argument(
        "--detailed-summary",
        "-D",
        action="store_true",
        help="Print a summary of all files created by the workflow. The "
        "has the following columns: filename, modification time, "
        "rule version, input file(s), shell command, status, plan.\n"
        "Thereby rule version contains the version "
        "the file was created with (see the version keyword of rules), and "
        "status denotes whether the file is missing, its input files are "
        "newer or if version or implementation of the rule changed since "
        "file creation. The input file and shell command columns are self "
        "explanatory. Finally the last column denotes whether the file "
        "will be updated or created during the next workflow execution.",
    )
    group_utils.add_argument(
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
        "formats are .tar, .tar.gz, .tar.bz2 and .tar.xz.",
    )
    group_utils.add_argument(
        "--cleanup-metadata",
        "--cm",
        nargs="+",
        metavar="FILE",
        help="Cleanup the metadata "
        "of given files. That means that snakemake removes any tracked "
        "version info, and any marks that files are incomplete.",
    )
    group_utils.add_argument(
        "--cleanup-shadow",
        action="store_true",
        help="Cleanup old shadow directories which have not been deleted due "
        "to failures or power loss.",
    )
    group_utils.add_argument(
        "--skip-script-cleanup",
        action="store_true",
        help="Don't delete wrapper scripts used for execution",
    )
    group_utils.add_argument(
        "--unlock", action="store_true", help="Remove a lock on the working directory."
    )
    group_utils.add_argument(
        "--list-version-changes",
        "--lv",
        action="store_true",
        help="List all output files that have been created with "
        "a different version (as determined by the version keyword).",
    )
    group_utils.add_argument(
        "--list-code-changes",
        "--lc",
        action="store_true",
        help="List all output files for which the rule body (run or shell) have "
        "changed in the Snakefile.",
    )
    group_utils.add_argument(
        "--list-input-changes",
        "--li",
        action="store_true",
        help="List all output files for which the defined input files have changed "
        "in the Snakefile (e.g. new input files were added in the rule "
        "definition or files were renamed). For listing input file "
        "modification in the filesystem, use --summary.",
    )
    group_utils.add_argument(
        "--list-params-changes",
        "--lp",
        action="store_true",
        help="List all output files for which the defined params have changed "
        "in the Snakefile.",
    )
    group_utils.add_argument(
        "--list-untracked",
        "--lu",
        action="store_true",
        help="List all files in the working directory that are not used in the  "
        "workflow. This can be used e.g. for identifying leftover files. Hidden files "
        "and directories are ignored.",
    )
    group_utils.add_argument(
        "--delete-all-output",
        action="store_true",
        help="Remove all files generated by the workflow. Use together with --dry-run "
        "to list files without actually deleting anything. Note that this will "
        "not recurse into subworkflows. Write-protected files are not removed. "
        "Nevertheless, use with care!",
    )
    group_utils.add_argument(
        "--delete-temp-output",
        action="store_true",
        help="Remove all temporary files generated by the workflow. Use together "
        "with --dry-run to list files without actually deleting anything. Note "
        "that this will not recurse into subworkflows.",
    )
    group_utils.add_argument(
        "--bash-completion",
        action="store_true",
        help="Output code to register bash completion for snakemake. Put the "
        "following in your .bashrc (including the accents): "
        "`snakemake --bash-completion` or issue it in an open terminal "
        "session.",
    )
    group_utils.add_argument(
        "--keep-incomplete",
        action="store_true",
        help="Do not remove incomplete output files by failed jobs.",
    )
    group_utils.add_argument(
        "--drop-metadata",
        action="store_true",
        help="Drop metadata file tracking information after job finishes. "
        "Provenance-information based reports (e.g. --report and the "
        "--list_x_changes functions) will be empty or incomplete.",
    )
    group_utils.add_argument("--version", "-v", action="version", version=__version__)

    group_output = parser.add_argument_group("OUTPUT")
    group_output.add_argument(
        "--reason",
        "-r",
        action="store_true",
        help="Print the reason for each executed rule.",
    )
    group_output.add_argument(
        "--gui",
        nargs="?",
        const="8000",
        metavar="PORT",
        type=str,
        help="Serve an HTML based user interface to the given network and "
        "port e.g. 168.129.10.15:8000. By default Snakemake is only "
        "available in the local network (default port: 8000). To make "
        "Snakemake listen to all ip addresses add the special host address "
        "0.0.0.0 to the url (0.0.0.0:8000). This is important if Snakemake "
        "is used in a virtualised environment like Docker. If possible, a "
        "browser window is opened.",
    )
    group_output.add_argument(
        "--printshellcmds",
        "-p",
        action="store_true",
        help="Print out the shell commands that will be executed.",
    )
    group_output.add_argument(
        "--debug-dag",
        action="store_true",
        help="Print candidate and selected jobs (including their wildcards) while "
        "inferring DAG. This can help to debug unexpected DAG topology or errors.",
    )
    group_output.add_argument(
        "--stats",
        metavar="FILE",
        help="Write stats about Snakefile execution in JSON format to the given file.",
    )
    group_output.add_argument(
        "--nocolor", action="store_true", help="Do not use a colored output."
    )
    group_output.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Do not output any progress or rule information.",
    )
    group_output.add_argument(
        "--print-compilation",
        action="store_true",
        help="Print the python representation of the workflow.",
    )

    group_output.add_argument(
        "--verbose", action="store_true", help="Print debugging output."
    )

    group_behavior = parser.add_argument_group("BEHAVIOR")
    group_behavior.add_argument(
        "--force-use-threads",
        dest="force_use_threads",
        action="store_true",
        help="Force threads rather than processes. Helpful if shared memory (/dev/shm) is full or unavailable.",
    )
    group_behavior.add_argument(
        "--allow-ambiguity",
        "-a",
        action="store_true",
        help=(
            "Don't check for ambiguous rules and simply use the first if "
            "several can produce the same file. This allows the user to "
            "prioritize rules by their order in the snakefile."
        ),
    )
    group_behavior.add_argument(
        "--nolock", action="store_true", help="Do not lock the working directory"
    )
    group_behavior.add_argument(
        "--ignore-incomplete",
        "--ii",
        action="store_true",
        help="Do not check for incomplete output files.",
    )
    group_behavior.add_argument(
        "--max-inventory-time",
        type=int,
        default=20,
        metavar="SECONDS",
        help="Spend at most SECONDS seconds to create a file inventory for the working directory. "
        "The inventory vastly speeds up file modification and existence checks when computing "
        "which jobs need to be executed. However, creating the inventory itself can be slow, e.g. on "
        "network file systems. Hence, we do not spend more than a given amount of time and fall back "
        "to individual checks for the rest.",
    )
    group_behavior.add_argument(
        "--latency-wait",
        "--output-wait",
        "-w",
        type=int,
        default=5,
        metavar="SECONDS",
        help="Wait given seconds if an output file of a job is not present after "
        "the job finished. This helps if your filesystem "
        "suffers from latency (default 5).",
    )
    group_behavior.add_argument(
        "--wait-for-files",
        nargs="*",
        metavar="FILE",
        help="Wait --latency-wait seconds for these "
        "files to be present before executing the workflow. "
        "This option is used internally to handle filesystem latency in cluster "
        "environments.",
    )
    group_behavior.add_argument(
        "--notemp",
        "--nt",
        action="store_true",
        help="Ignore temp() declarations. This is useful when running only "
        "a part of the workflow, since temp() would lead to deletion of "
        "probably needed files by other parts of the workflow.",
    )
    group_behavior.add_argument(
        "--keep-remote",
        action="store_true",
        help="Keep local copies of remote input files.",
    )
    group_behavior.add_argument(
        "--keep-target-files",
        action="store_true",
        help="Do not adjust the paths of given target files relative to the working directory.",
    )
    group_behavior.add_argument(
        "--allowed-rules",
        nargs="+",
        help="Only consider given rules. If omitted, all rules in Snakefile are "
        "used. Note that this is intended primarily for internal use and may "
        "lead to unexpected results otherwise.",
    )
    group_behavior.add_argument(
        "--max-jobs-per-second",
        default=10,
        type=float,
        help="Maximal number of cluster/drmaa jobs per second, default is 10, "
        "fractions allowed.",
    )
    group_behavior.add_argument(
        "--max-status-checks-per-second",
        default=10,
        type=float,
        help="Maximal number of job status checks per second, default is 10, "
        "fractions allowed.",
    )
    group_behavior.add_argument(
        "-T",
        "--restart-times",
        default=0,
        type=int,
        help="Number of times to restart failing jobs (defaults to 0).",
    )
    group_behavior.add_argument(
        "--attempt",
        default=1,
        type=int,
        help="Internal use only: define the initial value of the attempt "
        "parameter (default: 1).",
    )
    group_behavior.add_argument(
        "--wrapper-prefix",
        default="https://github.com/snakemake/snakemake-wrappers/raw/",
        help="Prefix for URL created from wrapper directive (default: "
        "https://github.com/snakemake/snakemake-wrappers/raw/). Set this to "
        "a different URL to use your fork or a local clone of the repository, "
        "e.g., use a git URL like 'git+file://path/to/your/local/clone@'.",
    )
    group_behavior.add_argument(
        "--default-remote-provider",
        choices=[
            "S3",
            "GS",
            "FTP",
            "SFTP",
            "S3Mocked",
            "gfal",
            "gridftp",
            "iRODS",
            "AzBlob",
        ],
        help="Specify default remote provider to be used for "
        "all input and output files that don't yet specify "
        "one.",
    )
    group_behavior.add_argument(
        "--default-remote-prefix",
        default="",
        help="Specify prefix for default remote provider. E.g. " "a bucket name.",
    )
    group_behavior.add_argument(
        "--no-shared-fs",
        action="store_true",
        help="Do not assume that jobs share a common file "
        "system. When this flag is activated, Snakemake will "
        "assume that the filesystem on a cluster node is not "
        "shared with other nodes. For example, this will lead "
        "to downloading remote files on each cluster node "
        "separately. Further, it won't take special measures "
        "to deal with filesystem latency issues. This option "
        "will in most cases only make sense in combination with "
        "--default-remote-provider. Further, when using --cluster "
        "you will have to also provide --cluster-status. "
        "Only activate this if you "
        "know what you are doing.",
    )
    group_behavior.add_argument(
        "--greediness",
        type=float,
        default=None,
        help="Set the greediness of scheduling. This value between 0 and 1 "
        "determines how careful jobs are selected for execution. The default "
        "value (1.0) provides the best speed and still acceptable scheduling "
        "quality.",
    )
    group_behavior.add_argument(
        "--no-hooks",
        action="store_true",
        help="Do not invoke onstart, onsuccess or onerror hooks after execution.",
    )
    group_behavior.add_argument(
        "--overwrite-shellcmd",
        help="Provide a shell command that shall be executed instead of those "
        "given in the workflow. "
        "This is for debugging purposes only.",
    )
    group_behavior.add_argument(
        "--debug",
        action="store_true",
        help="Allow to debug rules with e.g. PDB. This flag "
        "allows to set breakpoints in run blocks.",
    )
    group_behavior.add_argument(
        "--runtime-profile",
        metavar="FILE",
        help="Profile Snakemake and write the output to FILE. This requires yappi "
        "to be installed.",
    )
    group_behavior.add_argument(
        "--mode",
        choices=[Mode.default, Mode.subprocess, Mode.cluster],
        default=Mode.default,
        type=int,
        help="Set execution mode of Snakemake (internal use only).",
    )
    group_behavior.add_argument(
        "--show-failed-logs",
        action="store_true",
        help="Automatically display logs of failed jobs.",
    )
    group_behavior.add_argument(
        "--log-handler-script",
        metavar="FILE",
        default=None,
        help="Provide a custom script containing a function 'def log_handler(msg):'. "
        "Snakemake will call this function for every logging output (given as a dictionary msg)"
        "allowing to e.g. send notifications in the form of e.g. slack messages or emails.",
    )

    group_behavior.add_argument(
        "--log-service",
        default=None,
        choices=["none", "slack"],
        help="Set a specific messaging service for logging output."
        "Snakemake will notify the service on errors and completed execution."
        "Currently only slack is supported.",
    )

    group_cluster = parser.add_argument_group("CLUSTER")

    # TODO extend below description to explain the wildcards that can be used
    cluster_mode_group = group_cluster.add_mutually_exclusive_group()
    cluster_mode_group.add_argument(
        "--cluster",
        "-c",
        metavar="CMD",
        help=(
            "Execute snakemake rules with the given submit command, "
            "e.g. qsub. Snakemake compiles jobs into scripts that are "
            "submitted to the cluster with the given command, once all input "
            "files for a particular job are present.\n"
            "The submit command can be decorated to make it aware of certain "
            "job properties (name, rulename, input, output, params, wildcards, log, threads "
            "and dependencies (see the argument below)), e.g.:\n"
            "$ snakemake --cluster 'qsub -pe threaded {threads}'."
        ),
    ),
    cluster_mode_group.add_argument(
        "--cluster-sync",
        metavar="CMD",
        help=(
            "cluster submission command will block, returning the remote exit"
            "status upon remote termination (for example, this should be used"
            "if the cluster command is 'qsub -sync y' (SGE)"
        ),
    ),
    cluster_mode_group.add_argument(
        "--drmaa",
        nargs="?",
        const="",
        metavar="ARGS",
        help="Execute snakemake on a cluster accessed via DRMAA, "
        "Snakemake compiles jobs into scripts that are "
        "submitted to the cluster with the given command, once all input "
        "files for a particular job are present. ARGS can be used to "
        "specify options of the underlying cluster system, "
        "thereby using the job properties name, rulename, input, output, params, wildcards, log, "
        "threads and dependencies, e.g.: "
        "--drmaa ' -pe threaded {threads}'. Note that ARGS must be given in quotes and "
        "with a leading whitespace.",
    )

    group_cluster.add_argument(
        "--cluster-config",
        "-u",
        metavar="FILE",
        default=[],
        action="append",
        help=(
            "A JSON or YAML file that defines the wildcards used in 'cluster'"
            "for specific rules, instead of having them specified in the Snakefile. "
            "For example, for rule 'job' you may define: "
            "{ 'job' : { 'time' : '24:00:00' } } to specify the time for rule 'job'. "
            "You can specify more than one file.  The configuration files are merged "
            "with later values overriding earlier ones. This option is deprecated in favor "
            "of using --profile, see docs."
        ),
    ),
    group_cluster.add_argument(
        "--immediate-submit",
        "--is",
        action="store_true",
        help="Immediately submit all jobs to the cluster instead of waiting "
        "for present input files. This will fail, unless you make "
        "the cluster aware of job dependencies, e.g. via:\n"
        "$ snakemake --cluster 'sbatch --dependency {dependencies}.\n"
        "Assuming that your submit script (here sbatch) outputs the "
        "generated job id to the first stdout line, {dependencies} will "
        "be filled with space separated job ids this job depends on.",
    )
    group_cluster.add_argument(
        "--jobscript",
        "--js",
        metavar="SCRIPT",
        help="Provide a custom job script for submission to the cluster. "
        "The default script resides as 'jobscript.sh' in the "
        "installation directory.",
    )
    group_cluster.add_argument(
        "--jobname",
        "--jn",
        default="snakejob.{name}.{jobid}.sh",
        metavar="NAME",
        help="Provide a custom name for the jobscript that is submitted to the "
        'cluster (see --cluster). NAME is "snakejob.{name}.{jobid}.sh" '
        "per default. The wildcard {jobid} has to be present in the name.",
    )
    group_cluster.add_argument(
        "--cluster-status",
        help="Status command for cluster execution. This is only considered "
        "in combination with the --cluster flag. If provided, Snakemake will "
        "use the status command to determine if a job has finished successfully "
        "or failed. For this it is necessary that the submit command provided "
        "to --cluster returns the cluster job id. Then, the status command "
        "will be invoked with the job id. Snakemake expects it to return "
        "'success' if the job was successfull, 'failed' if the job failed and "
        "'running' if the job still runs.",
    )
    group_cluster.add_argument(
        "--drmaa-log-dir",
        metavar="DIR",
        help="Specify a directory in which stdout and stderr files of DRMAA"
        " jobs will be written. The value may be given as a relative path,"
        " in which case Snakemake will use the current invocation directory"
        " as the origin. If given, this will override any given '-o' and/or"
        " '-e' native specification. If not given, all DRMAA stdout and"
        " stderr files are written to the current working directory.",
    )

    group_cloud = parser.add_argument_group("CLOUD")
    group_kubernetes = parser.add_argument_group("KUBERNETES")
    group_tibanna = parser.add_argument_group("TIBANNA")
    group_google_life_science = parser.add_argument_group("GOOGLE_LIFE_SCIENCE")
    group_tes = parser.add_argument_group("TES")

    group_kubernetes.add_argument(
        "--kubernetes",
        metavar="NAMESPACE",
        nargs="?",
        const="default",
        help="Execute workflow in a kubernetes cluster (in the cloud). "
        "NAMESPACE is the namespace you want to use for your job (if nothing "
        "specified: 'default'). "
        "Usually, this requires --default-remote-provider and "
        "--default-remote-prefix to be set to a S3 or GS bucket where your . "
        "data shall be stored. It is further advisable to activate conda "
        "integration via --use-conda.",
    )
    group_kubernetes.add_argument(
        "--container-image",
        metavar="IMAGE",
        help="Docker image to use, e.g., when submitting jobs to kubernetes "
        "Defaults to 'https://hub.docker.com/r/snakemake/snakemake', tagged with "
        "the same version as the currently running Snakemake instance. "
        "Note that overwriting this value is up to your responsibility. "
        "Any used image has to contain a working snakemake installation "
        "that is compatible with (or ideally the same as) the currently "
        "running version.",
    )

    group_tibanna.add_argument(
        "--tibanna",
        action="store_true",
        help="Execute workflow on AWS cloud using Tibanna. This requires "
        "--default-remote-prefix to be set to S3 bucket name and prefix"
        " (e.g. 'bucketname/subdirectory') where input is already stored"
        " and output will be sent to. Using --tibanna implies --default-resources"
        " is set as default. Optionally, use --precommand to"
        " specify any preparation command to run before snakemake command"
        " on the cloud (inside snakemake container on Tibanna VM)."
        " Also, --use-conda, --use-singularity, --config, --configfile are"
        " supported and will be carried over.",
    )
    group_tibanna.add_argument(
        "--tibanna-sfn",
        help="Name of Tibanna Unicorn step function (e.g. tibanna_unicorn_monty)."
        "This works as serverless scheduler/resource allocator and must be "
        "deployed first using tibanna cli. (e.g. tibanna deploy_unicorn --usergroup="
        "monty --buckets=bucketname)",
    )
    group_tibanna.add_argument(
        "--precommand",
        help="Any command to execute before snakemake command on AWS cloud "
        "such as wget, git clone, unzip, etc. This is used with --tibanna."
        "Do not include input/output download/upload commands - file transfer"
        " between S3 bucket and the run environment (container) is automatically"
        " handled by Tibanna.",
    )
    group_tibanna.add_argument(
        "--tibanna-config",
        nargs="+",
        help="Additional tibanna config e.g. --tibanna-config spot_instance=true subnet="
        "<subnet_id> security group=<security_group_id>",
    )
    group_google_life_science.add_argument(
        "--google-lifesciences",
        action="store_true",
        help="Execute workflow on Google Cloud cloud using the Google Life. "
        " Science API. This requires default application credentials (json) "
        " to be created and export to the environment to use Google Cloud "
        " Storage, Compute Engine, and Life Sciences. The credential file "
        " should be exported as GOOGLE_APPLICATION_CREDENTIALS for snakemake "
        " to discover. Also, --use-conda, --use-singularity, --config, "
        "--configfile are supported and will be carried over.",
    )
    group_google_life_science.add_argument(
        "--google-lifesciences-regions",
        nargs="+",
        default=["us-east1", "us-west1", "us-central1"],
        help="Specify one or more valid instance regions (defaults to US)",
    )
    group_google_life_science.add_argument(
        "--google-lifesciences-location",
        help="The Life Sciences API service used to schedule the jobs. "
        " E.g., us-centra1 (Iowa) and europe-west2 (London) "
        " Watch the terminal output to see all options found to be available. "
        " If not specified, defaults to the first found with a matching prefix "
        " from regions specified with --google-lifesciences-regions.",
    )
    group_google_life_science.add_argument(
        "--google-lifesciences-keep-cache",
        action="store_true",
        help="Cache workflows in your Google Cloud Storage Bucket specified "
        "by --default-remote-prefix/{source}/{cache}. Each workflow working "
        "directory is compressed to a .tar.gz, named by the hash of the "
        "contents, and kept in Google Cloud Storage. By default, the caches "
        "are deleted at the shutdown step of the workflow.",
    )

    group_tes.add_argument(
        "--tes",
        metavar="URL",
        help="Send workflow tasks to GA4GH TES server specified by url.",
    )

    group_conda = parser.add_argument_group("CONDA")

    group_conda.add_argument(
        "--use-conda",
        action="store_true",
        help="If defined in the rule, run job in a conda environment. "
        "If this flag is not set, the conda directive is ignored.",
    )
    group_conda.add_argument(
        "--list-conda-envs",
        action="store_true",
        help="List all conda environments and their location on " "disk.",
    )
    group_conda.add_argument(
        "--conda-prefix",
        metavar="DIR",
        help="Specify a directory in which the 'conda' and 'conda-archive' "
        "directories are created. These are used to store conda environments "
        "and their archives, respectively. If not supplied, the value is set "
        "to the '.snakemake' directory relative to the invocation directory. "
        "If supplied, the `--use-conda` flag must also be set. The value may "
        "be given as a relative path, which will be extrapolated to the "
        "invocation directory, or as an absolute path.",
    )
    group_conda.add_argument(
        "--conda-cleanup-envs",
        action="store_true",
        help="Cleanup unused conda environments.",
    )

    from snakemake.deployment.conda import CondaCleanupMode

    group_conda.add_argument(
        "--conda-cleanup-pkgs",
        type=CondaCleanupMode,
        const=CondaCleanupMode.tarballs,
        choices=list(CondaCleanupMode),
        nargs="?",
        help="Cleanup conda packages after creating environments. "
        "In case of 'tarballs' mode, will clean up all downloaded package tarballs. "
        "In case of 'cache' mode, will additionally clean up unused package caches. "
        "If mode is omitted, will default to only cleaning up the tarballs.",
    )
    group_conda.add_argument(
        "--conda-create-envs-only",
        action="store_true",
        help="If specified, only creates the job-specific "
        "conda environments then exits. The `--use-conda` "
        "flag must also be set.",
    )
    group_conda.add_argument(
        "--conda-frontend",
        default="conda",
        choices=["conda", "mamba"],
        help="Choose the conda frontend for installing environments. "
        "Caution: mamba is much faster, but still in beta test.",
    )

    group_singularity = parser.add_argument_group("SINGULARITY")

    group_singularity.add_argument(
        "--use-singularity",
        action="store_true",
        help="If defined in the rule, run job within a singularity container. "
        "If this flag is not set, the singularity directive is ignored.",
    )
    group_singularity.add_argument(
        "--singularity-prefix",
        metavar="DIR",
        help="Specify a directory in which singularity images will be stored."
        "If not supplied, the value is set "
        "to the '.snakemake' directory relative to the invocation directory. "
        "If supplied, the `--use-singularity` flag must also be set. The value "
        "may be given as a relative path, which will be extrapolated to the "
        "invocation directory, or as an absolute path.",
    )
    group_singularity.add_argument(
        "--singularity-args",
        default="",
        metavar="ARGS",
        help="Pass additional args to singularity.",
    )

    group_env_modules = parser.add_argument_group("ENVIRONMENT MODULES")

    group_env_modules.add_argument(
        "--use-envmodules",
        action="store_true",
        help="If defined in the rule, run job within the given environment "
        "modules, loaded in the given order. This can be combined with "
        "--use-conda and --use-singularity, which will then be only used as a "
        "fallback for rules which don't define environment modules.",
    )

    return parser


def main(argv=None):
    """Main entry point."""

    if sys.version_info < MIN_PY_VERSION:
        print(
            "Snakemake requires at least Python {}.".format(MIN_PY_VERSION),
            file=sys.stderr,
        )
        exit(1)

    parser = get_argument_parser()
    args = parser.parse_args(argv)

    if args.profile:
        # reparse args while inferring config file from profile
        parser = get_argument_parser(args.profile)
        args = parser.parse_args(argv)

        def adjust_path(f):
            if os.path.exists(f) or os.path.isabs(f):
                return f
            else:
                return get_profile_file(args.profile, f, return_default=True)

        # update file paths to be relative to the profile
        # (if they do not exist relative to CWD)
        if args.jobscript:
            args.jobscript = adjust_path(args.jobscript)
        if args.cluster:
            args.cluster = adjust_path(args.cluster)
        if args.cluster_sync:
            args.cluster_sync = adjust_path(args.cluster_sync)
        if args.cluster_status:
            args.cluster_status = adjust_path(args.cluster_status)
        if args.report_stylesheet:
            args.report_stylesheet = adjust_path(args.report_stylesheet)

    if args.bash_completion:
        cmd = b"complete -o bashdefault -C snakemake-bash-completion snakemake"
        sys.stdout.buffer.write(cmd)
        sys.exit(0)

    if args.batch is not None and args.forceall:
        print(
            "--batch may not be combined with --forceall, because recomputed upstream "
            "jobs in subsequent batches may render already obtained results outdated."
        )

    try:
        resources = parse_resources(args.resources)
        config = parse_config(args)

        # Cloud executors should have default-resources flag
        if (
            (args.default_resources is not None and not args.default_resources)
            or (args.tibanna and not args.default_resources)
            or (args.google_lifesciences and not args.default_resources)
        ):
            args.default_resources = [
                "mem_mb=max(2*input.size_mb, 1000)",
                "disk_mb=max(2*input.size_mb, 1000)",
            ]
        default_resources = DefaultResources(args.default_resources)
        batch = parse_batch(args)
        overwrite_threads = parse_set_threads(args)

        overwrite_scatter = parse_set_scatter(args)

        overwrite_groups = parse_groups(args)
        group_components = parse_group_components(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        print("", file=sys.stderr)
        sys.exit(1)

    local_exec = not (
        args.print_compilation
        or args.cluster
        or args.cluster_sync
        or args.drmaa
        or args.google_lifesciences
        or args.kubernetes
        or args.tibanna
        or args.tes
        or args.list_code_changes
        or args.list_conda_envs
        or args.list_input_changes
        or args.list_params_changes
        or args.list
        or args.list_target_rules
        or args.list_untracked
        or args.list_version_changes
        or args.export_cwl
        or args.generate_unit_tests
        or args.dag
        or args.d3dag
        or args.filegraph
        or args.rulegraph
        or args.summary
        or args.lint
        or args.report
        or args.gui
        or args.archive
    )

    if args.cores is not None:
        if args.cores == "all":
            args.cores = available_cpu_count()
        else:
            try:
                args.cores = int(args.cores)
            except ValueError:
                print(
                    "Error parsing number of cores (--cores, --jobs, -j): must be integer, empty, or 'all'.",
                    file=sys.stderr,
                )
                sys.exit(1)
    if args.cluster or args.cluster_sync or args.drmaa:
        if args.cores is None:
            if args.dryrun:
                args.cores = 1
            else:
                print(
                    "Error: you need to specify the maximum number of jobs to "
                    "be queued or executed at the same time with --jobs.",
                    file=sys.stderr,
                )
                sys.exit(1)
    elif args.cores is None:
        if local_exec and not (args.dryrun or args.unlock):
            print(
                "Error: you need to specify the maximum number of CPU cores to "
                "be used at the same time. If you want to use N cores, say --cores N or "
                "-jN. For all cores on your system (be sure that this is appropriate) "
                "use --cores all. For no parallelization use --cores 1 or -j1.",
                file=sys.stderr,
            )
            sys.exit(1)
        else:
            args.cores = 1

    if args.drmaa_log_dir is not None:
        if not os.path.isabs(args.drmaa_log_dir):
            args.drmaa_log_dir = os.path.abspath(os.path.expanduser(args.drmaa_log_dir))

    if args.runtime_profile:
        import yappi

        yappi.start()

    if args.immediate_submit and not args.notemp:
        print(
            "Error: --immediate-submit has to be combined with --notemp, "
            "because temp file handling is not supported in this mode.",
            file=sys.stderr,
        )
        sys.exit(1)

    if (args.conda_prefix or args.conda_create_envs_only) and not args.use_conda:
        print(
            "Error: --use-conda must be set if --conda-prefix or "
            "--create-envs-only is set.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.singularity_prefix and not args.use_singularity:
        print(
            "Error: --use_singularity must be set if --singularity-prefix " "is set.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.kubernetes and (
        not args.default_remote_provider or not args.default_remote_prefix
    ):
        print(
            "Error: --kubernetes must be combined with "
            "--default-remote-provider and --default-remote-prefix, see "
            "https://snakemake.readthedocs.io/en/stable/executing/cloud.html"
            "#executing-a-snakemake-workflow-via-kubernetes",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.tibanna:
        if not args.default_remote_prefix:
            print(
                "Error: --tibanna must be combined with --default-remote-prefix "
                "to provide bucket name and subdirectory (prefix) "
                "(e.g. 'bucketname/projectname'",
                file=sys.stderr,
            )
            sys.exit(1)
        args.default_remote_prefix = args.default_remote_prefix.rstrip("/")
        if not args.tibanna_sfn:
            args.tibanna_sfn = os.environ.get("TIBANNA_DEFAULT_STEP_FUNCTION_NAME", "")
            if not args.tibanna_sfn:
                print(
                    "Error: to use --tibanna, either --tibanna-sfn or environment variable "
                    "TIBANNA_DEFAULT_STEP_FUNCTION_NAME must be set and exported "
                    "to provide name of the tibanna unicorn step function "
                    "(e.g. 'tibanna_unicorn_monty'). The step function must be deployed first "
                    "using tibanna cli (e.g. tibanna deploy_unicorn --usergroup=monty "
                    "--buckets=bucketname)",
                    file=sys.stderr,
                )
                sys.exit(1)

    if args.google_lifesciences:
        if not os.environ.get("GOOGLE_APPLICATION_CREDENTIALS"):
            print(
                "Error: GOOGLE_APPLICATION_CREDENTIALS environment variable must "
                "be available for --google-lifesciences",
                file=sys.stderr,
            )
            sys.exit(1)

        if not args.default_remote_prefix:
            print(
                "Error: --google-life-sciences must be combined with "
                " --default-remote-prefix to provide bucket name and "
                "subdirectory (prefix) (e.g. 'bucketname/projectname'",
                file=sys.stderr,
            )
            sys.exit(1)

    if args.delete_all_output and args.delete_temp_output:
        print(
            "Error: --delete-all-output and --delete-temp-output are mutually exclusive.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.snakefile is None:
        for p in SNAKEFILE_CHOICES:
            if os.path.exists(p):
                args.snakefile = p
                break
        if args.snakefile is None:
            print(
                "Error: no Snakefile found, tried {}.".format(
                    ", ".join(SNAKEFILE_CHOICES), file=sys.stderr
                )
            )
            sys.exit(1)

    if args.gui is not None:
        try:
            import snakemake.gui as gui
        except ImportError:
            print(
                "Error: GUI needs Flask to be installed. Install "
                "with easy_install or contact your administrator.",
                file=sys.stderr,
            )
            sys.exit(1)

        _logging.getLogger("werkzeug").setLevel(_logging.ERROR)

        _snakemake = partial(snakemake, os.path.abspath(args.snakefile))
        gui.register(_snakemake, args)

        if ":" in args.gui:
            host, port = args.gui.split(":")
        else:
            port = args.gui
            host = "127.0.0.1"

        url = "http://{}:{}".format(host, port)
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
            gui.app.run(debug=False, threaded=True, port=int(port), host=host)

        except (KeyboardInterrupt, SystemExit):
            # silently close
            pass
    else:
        log_handler = []
        if args.log_handler_script is not None:
            if not os.path.exists(args.log_handler_script):
                print(
                    "Error: no log handler script found, {}.".format(
                        args.log_handler_script
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)
            log_script = SourceFileLoader("log", args.log_handler_script).load_module()
            try:
                log_handler.append(log_script.log_handler)
            except:
                print(
                    'Error: Invalid log handler script, {}. Expect python function "log_handler(msg)".'.format(
                        args.log_handler_script
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)

        if args.log_service == "slack":
            slack_logger = logging.SlackLogger()
            log_handler.append(slack_logger.log_handler)

        if args.edit_notebook:
            from snakemake import notebook

            args.target = [args.edit_notebook]
            args.force = True
            args.edit_notebook = notebook.Listen(args.notebook_listen)

        success = snakemake(
            args.snakefile,
            batch=batch,
            cache=args.cache,
            report=args.report,
            report_stylesheet=args.report_stylesheet,
            lint=args.lint,
            generate_unit_tests=args.generate_unit_tests,
            listrules=args.list,
            list_target_rules=args.list_target_rules,
            cores=args.cores,
            local_cores=args.local_cores,
            nodes=args.cores,
            resources=resources,
            overwrite_threads=overwrite_threads,
            overwrite_scatter=overwrite_scatter,
            default_resources=default_resources,
            config=config,
            configfiles=args.configfile,
            config_args=args.config,
            workdir=args.directory,
            targets=args.target,
            dryrun=args.dryrun,
            printshellcmds=args.printshellcmds,
            printreason=args.reason,
            debug_dag=args.debug_dag,
            printdag=args.dag,
            printrulegraph=args.rulegraph,
            printfilegraph=args.filegraph,
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
            kubernetes=args.kubernetes,
            container_image=args.container_image,
            tibanna=args.tibanna,
            tibanna_sfn=args.tibanna_sfn,
            google_lifesciences=args.google_lifesciences,
            google_lifesciences_regions=args.google_lifesciences_regions,
            google_lifesciences_location=args.google_lifesciences_location,
            google_lifesciences_cache=args.google_lifesciences_keep_cache,
            tes=args.tes,
            precommand=args.precommand,
            preemption_default=args.preemption_default,
            preemptible_rules=args.preemptible_rules,
            tibanna_config=args.tibanna_config,
            jobname=args.jobname,
            immediate_submit=args.immediate_submit,
            standalone=True,
            ignore_ambiguity=args.allow_ambiguity,
            lock=not args.nolock,
            unlock=args.unlock,
            cleanup_metadata=args.cleanup_metadata,
            conda_cleanup_envs=args.conda_cleanup_envs,
            cleanup_shadow=args.cleanup_shadow,
            cleanup_scripts=not args.skip_script_cleanup,
            force_incomplete=args.rerun_incomplete,
            ignore_incomplete=args.ignore_incomplete,
            list_version_changes=args.list_version_changes,
            list_code_changes=args.list_code_changes,
            list_input_changes=args.list_input_changes,
            list_params_changes=args.list_params_changes,
            list_untracked=args.list_untracked,
            summary=args.summary,
            detailed_summary=args.detailed_summary,
            archive=args.archive,
            delete_all_output=args.delete_all_output,
            delete_temp_output=args.delete_temp_output,
            print_compilation=args.print_compilation,
            verbose=args.verbose,
            debug=args.debug,
            jobscript=args.jobscript,
            notemp=args.notemp,
            keep_remote_local=args.keep_remote,
            greediness=args.greediness,
            no_hooks=args.no_hooks,
            overwrite_shellcmd=args.overwrite_shellcmd,
            latency_wait=args.latency_wait,
            wait_for_files=args.wait_for_files,
            keep_target_files=args.keep_target_files,
            allowed_rules=args.allowed_rules,
            max_jobs_per_second=args.max_jobs_per_second,
            max_status_checks_per_second=args.max_status_checks_per_second,
            restart_times=args.restart_times,
            attempt=args.attempt,
            force_use_threads=args.force_use_threads,
            use_conda=args.use_conda,
            conda_frontend=args.conda_frontend,
            conda_prefix=args.conda_prefix,
            conda_cleanup_pkgs=args.conda_cleanup_pkgs,
            list_conda_envs=args.list_conda_envs,
            use_singularity=args.use_singularity,
            use_env_modules=args.use_envmodules,
            singularity_prefix=args.singularity_prefix,
            shadow_prefix=args.shadow_prefix,
            singularity_args=args.singularity_args,
            scheduler=args.scheduler,
            scheduler_ilp_solver=args.scheduler_ilp_solver,
            conda_create_envs_only=args.conda_create_envs_only,
            mode=args.mode,
            wrapper_prefix=args.wrapper_prefix,
            default_remote_provider=args.default_remote_provider,
            default_remote_prefix=args.default_remote_prefix,
            assume_shared_fs=not args.no_shared_fs,
            cluster_status=args.cluster_status,
            export_cwl=args.export_cwl,
            show_failed_logs=args.show_failed_logs,
            wms_monitor=args.wms_monitor,
            keep_incomplete=args.keep_incomplete,
            keep_metadata=not args.drop_metadata,
            edit_notebook=args.edit_notebook,
            envvars=args.envvars,
            overwrite_groups=overwrite_groups,
            group_components=group_components,
            max_inventory_wait_time=args.max_inventory_time,
            log_handler=log_handler,
            execute_subworkflows=not args.no_subworkflows,
        )

    if args.runtime_profile:
        with open(args.runtime_profile, "w") as out:
            profile = yappi.get_func_stats()
            profile.sort("totaltime")
            profile.print_all(
                out=out,
                columns={
                    0: ("name", 120),
                    1: ("ncall", 10),
                    2: ("tsub", 8),
                    3: ("ttot", 8),
                    4: ("tavg", 8),
                },
            )

    sys.exit(0 if success else 1)


def bash_completion(snakefile="Snakefile"):
    """Entry point for bash completion."""
    if not len(sys.argv) >= 2:
        print(
            "Calculate bash completion for snakemake. This tool shall not be invoked by hand."
        )
        sys.exit(1)

    def print_candidates(candidates):
        if candidates:
            candidates = sorted(set(candidates))
            ## Use bytes for avoiding '^M' under Windows.
            sys.stdout.buffer.write(b"\n".join(s.encode() for s in candidates))

    prefix = sys.argv[2]

    if prefix.startswith("-"):
        print_candidates(
            action.option_strings[0]
            for action in get_argument_parser()._actions
            if action.option_strings and action.option_strings[0].startswith(prefix)
        )
    else:
        candidates = []
        files = glob.glob("{}*".format(prefix))
        if files:
            candidates.extend(files)
        if os.path.exists(snakefile):
            workflow = Workflow(snakefile=snakefile)
            workflow.include(snakefile)

            candidates.extend(
                [file for file in workflow.concrete_files if file.startswith(prefix)]
                + [rule.name for rule in workflow.rules if rule.name.startswith(prefix)]
            )
        if len(candidates) > 0:
            print_candidates(candidates)
    sys.exit(0)
