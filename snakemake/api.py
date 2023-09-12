__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import sys

from snakemake.common import MIN_PY_VERSION

if sys.version_info < MIN_PY_VERSION:
    raise ValueError(f"Snakemake requires at least Python {'.'.join(MIN_PY_VERSION)}.")

import os
from functools import partial
import importlib

from snakemake_interface_executor_plugins.utils import ExecMode

from snakemake.workflow import Workflow
from snakemake.exceptions import (
    print_exception,
    WorkflowError,
)
from snakemake.logging import setup_logger, logger
from snakemake.io import load_configfile
from snakemake.shell import shell
from snakemake.utils import update_config
from snakemake.common import (
    MIN_PY_VERSION,
    RERUN_TRIGGERS,
    __version__,
    dict_to_key_value_args,
)
from snakemake.resources import DefaultResources


def snakemake(
    snakefile,
    batch=None,
    cache=None,
    report=None,
    report_stylesheet=None,
    containerize=False,
    lint=None,
    generate_unit_tests=None,
    listrules=False,
    list_target_rules=False,
    cores=1,
    nodes=None,
    local_cores=1,
    max_threads=None,
    resources=dict(),
    overwrite_threads=None,
    overwrite_scatter=None,
    overwrite_resource_scopes=None,
    default_resources=None,
    overwrite_resources=None,
    config=dict(),
    configfiles=None,
    config_args=None,
    workdir=None,
    targets=None,
    target_jobs=None,
    dryrun=False,
    touch=False,
    forcetargets=False,
    forceall=False,
    forcerun=[],
    until=[],
    omit_from=[],
    prioritytargets=[],
    stats=None,
    printshellcmds=False,
    debug_dag=False,
    printdag=False,
    printrulegraph=False,
    printfilegraph=False,
    printd3dag=False,
    nocolor=False,
    quiet=False,
    keepgoing=False,
    slurm=None,
    slurm_jobstep=None,
    rerun_triggers=RERUN_TRIGGERS,
    cluster=None,
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
    cleanup_containers=False,
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
    all_temp=False,
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
    mode=ExecMode.default,
    wrapper_prefix=None,
    kubernetes=None,
    container_image=None,
    k8s_cpu_scalar=1.0,
    k8s_service_account_name=None,
    flux=False,
    tibanna=False,
    tibanna_sfn=None,
    az_batch=False,
    az_batch_enable_autoscale=False,
    az_batch_account_url=None,
    google_lifesciences=False,
    google_lifesciences_regions=None,
    google_lifesciences_location=None,
    google_lifesciences_cache=False,
    google_lifesciences_service_account_email=None,
    google_lifesciences_network=None,
    google_lifesciences_subnetwork=None,
    tes=None,
    preemption_default=None,
    preemptible_rules=None,
    precommand="",
    default_remote_provider=None,
    default_remote_prefix="",
    tibanna_config=False,
    assume_shared_fs=True,
    cluster_status=None,
    cluster_cancel=None,
    cluster_cancel_nargs=None,
    cluster_sidecar=None,
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
    conda_not_block_search_path_envvars=False,
    scheduler_solver_path=None,
    conda_base_path=None,
    local_groupid="local",
    executor_args=None,
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
        workdir (str):              path to the working directory (default None)
        targets (list):             list of targets, e.g. rule or file names (default None)
        target_jobs (dict):         list of snakemake.target_jobs.TargetSpec objects directly targeting specific jobs (default None)
        dryrun (bool):              only dry-run the workflow (default False)
        touch (bool):               only touch all output files if present (default False)
        forcetargets (bool):        force given targets to be re-created (default False)
        forceall (bool):            force all output files to be re-created (default False)
        forcerun (list):            list of files and rules that shall be re-created/re-executed (default [])
        execute_subworkflows (bool):   execute subworkflows if present (default True)
        prioritytargets (list):     list of targets that shall be run with maximum priority (default [])
        stats (str):                path to file that shall contain stats about the workflow execution (default None)
        printshellcmds (bool):      print the shell command of each job (default False)
        printdag (bool):            print the dag in the graphviz dot language (default False)
        printrulegraph (bool):      print the graph of rules in the graphviz dot language (default False)
        printfilegraph (bool):      print the graph of rules with their input and output files in the graphviz dot language (default False)
        printd3dag (bool):          print a D3.js compatible JSON representation of the DAG (default False)
        nocolor (bool):             do not print colored output (default False)
        quiet (bool):               do not print any default job information (default False)
        keepgoing (bool):           keep going upon errors (default False)
        cluster (str):              submission command of a cluster or batch system to use, e.g. qsub (default None)
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
        cleanup_containers (bool):  delete unused (singularity) containers (default False)
        force_incomplete (bool):    force the re-creation of incomplete files (default False)
        ignore_incomplete (bool):   ignore incomplete files (default False)
        list_version_changes (bool): list output files with changed rule version (default False)
        list_code_changes (bool):   list output files with changed rule code (default False)
        list_input_changes (bool):  list output files with changed input files (default False)
        list_params_changes (bool): list output files with changed params (default False)
        list_untracked (bool):      list files in the workdir that are not used in the workflow (default False)
        summary (bool):             list summary of all output files and their status (default False)
        archive (str):              archive workflow into the given tarball
        delete_all_output (bool):    remove all files generated by the workflow (default False)
        delete_temp_output (bool):   remove all temporary files generated by the workflow (default False)
        latency_wait (int):         how many seconds to wait for an output file to appear after the execution of a job, e.g. to handle filesystem latency (default 3)
        wait_for_files (list):      wait for given files to be present before executing the workflow
        list_resources (bool):      list resources used in the workflow (default False)
        summary (bool):             list summary of all output files and their status (default False). If no option is specified a basic summary will be output. If 'detailed' is added as an option e.g --summary detailed, extra info about the input and shell commands will be included
        detailed_summary (bool):    list summary of all input and output files and their status (default False)
        print_compilation (bool):   print the compilation of the snakefile (default False)
        debug (bool):               allow to use the debugger within rules
        notemp (bool):              ignore temp file flags, e.g. do not delete output files marked as a temp after use (default False)
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
        force_use_threads:          whether to force the use of threads over processes. helpful if shared memory is full or unavailable (default False)
        use_conda (bool):           use conda environments for each job (defined with conda directive of rules)
        use_singularity (bool):     run jobs in singularity containers (if defined with singularity directive)
        use_env_modules (bool):     load environment modules if defined in rules
        singularity_args (str):     additional arguments to pass to a singularity
        conda_prefix (str):         the directory in which conda environments will be created (default None)
        conda_cleanup_pkgs (snakemake.deployment.conda.CondaCleanupMode):
                                    whether to clean up conda tarballs after env creation (default None), valid values: "tarballs", "cache"
        singularity_prefix (str):   the directory to which singularity images will be pulled (default None)
        shadow_prefix (str):        prefix for shadow directories. The job-specific shadow directories will be created in $SHADOW_PREFIX/shadow/ (default None)
        conda_create_envs_only (bool):    if specified, only builds the conda environments specified for each job, then exits.
        list_conda_envs (bool):     list conda environments and their location on disk.
        mode (snakemake.common.Mode): execution mode
        wrapper_prefix (str):       prefix for wrapper script URLs (default None)
        kubernetes (str):           submit jobs to Kubernetes, using the given namespace.
        container_image (str):      Docker image to use, e.g., for Kubernetes.
        k8s_cpu_scalar (float):     What proportion of each k8s node's CPUs are availabe to snakemake?
        k8s_service_account_name (str): Custom k8s service account, needed for workload identity.
        flux (bool):                Launch workflow to flux cluster.
        default_remote_provider (str): default remote provider to use instead of local files (e.g. S3, GS)
        default_remote_prefix (str): prefix for default remote provider (e.g. name of the bucket).
        tibanna (bool):             submit jobs to AWS cloud using Tibanna.
        tibanna_sfn (str):          Step function (Unicorn) name of Tibanna (e.g. tibanna_unicorn_monty). This must be deployed first using tibanna cli.
        az_batch (bool):            Submit jobs to azure batch.
        az_batch_enable_autoscale (bool): Enable autoscaling of the azure batch pool nodes. This sets the initial dedicated node pool count to zero and resizes the pool only after 5 minutes. So this flag is only recommended for relatively long running jobs.,
        az_batch_account_url (str): Azure batch account url.
        google_lifesciences (bool): submit jobs to Google Cloud Life Sciences (pipelines API).
        google_lifesciences_regions (list): a list of regions (e.g., us-east1)
        google_lifesciences_location (str): Life Sciences API location (e.g., us-central1)
        google_lifesciences_cache (bool): save a cache of the compressed working directories in Google Cloud Storage for later usage.
        google_lifesciences_service_account_email (str): Service account to install on Google pipelines API VM instance.
        google_lifesciences_network (str): Network name for Google VM instances.
        google_lifesciences_subnetwork (str): Subnetwork name for Google VM instances.
        tes (str):                  Execute workflow tasks on GA4GH TES server given by URL.
        precommand (str):           commands to run on AWS cloud before the snakemake command (e.g. wget, git clone, unzip, etc). Use with --tibanna.
        preemption_default (int):   set a default number of preemptible instance retries (for Google Life Sciences executor only)
        preemptible_rules (list):    define custom preemptible instance retries for specific rules (for Google Life Sciences executor only)
        tibanna_config (list):      Additional tibanna config e.g. --tibanna-config spot_instance=true subnet=<subnet_id> security group=<security_group_id>
        assume_shared_fs (bool):    assume that cluster nodes share a common filesystem (default true).
        cluster_status (str):       status command for cluster execution. If None, Snakemake will rely on flag files. Otherwise, it expects the command to return "success", "failure" or "running" when executing with a cluster jobid as a single argument.
        cluster_cancel (str):       command to cancel multiple job IDs (like SLURM 'scancel') (default None)
        cluster_cancel_nargs (int): maximal number of job ids to pass to cluster_cancel (default 1000)
        cluster_sidecar (str):      command that starts a sidecar process, see cluster documentation (default None)
        export_cwl (str):           Compile workflow to CWL and save to given file
        log_handler (function):     redirect snakemake output to this custom log handler, a function that takes a log message dictionary (see below) as its only argument (default None). The log message dictionary for the log handler has to following entries:
        keep_incomplete (bool):     keep incomplete output files of failed jobs
        edit_notebook (object):     "notebook.EditMode" object to configure notebook server for interactive editing of a rule notebook. If None, do not edit.
        scheduler (str):            Select scheduling algorithm (default ilp)
        scheduler_ilp_solver (str): Set solver for ilp scheduler.
        overwrite_groups (dict):    Rule to group assignments (default None)
        group_components (dict):    Number of connected components given groups shall span before being split up (1 by default if empty)
        conda_not_block_search_path_envvars (bool): Do not block search path envvars (R_LIBS, PYTHONPATH, ...) when using conda environments.
        scheduler_solver_path (str): Path to Snakemake environment (this can be used to e.g. overwrite the search path for the ILP solver used during scheduling).
        conda_base_path (str):      Path to conda base environment (this can be used to overwrite the search path for conda, mamba, and activate).
        local_groupid (str):        Local groupid to use as a placeholder for groupid-referrring input functions of local jobs (internal use only, default: local).
        log_handler (list):         redirect snakemake output to this list of custom log handlers, each a function that takes a log message dictionary (see below) as its only argument (default []). The log message dictionary for the log handler has to following entries:
        executor_args (dataclasses.Dataclass):  custom Data class to pass to custom executors for more flexibility
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

    # Azure batch uses compute engine and storage
    if az_batch:
        assume_shared_fs = False
        default_remote_provider = "AzBlob"

    # Google Cloud Life Sciences API uses compute engine and storage
    if google_lifesciences:
        assume_shared_fs = False
        default_remote_provider = "GS"
        default_remote_prefix = default_remote_prefix.rstrip("/")
    if kubernetes:
        assume_shared_fs = False

    # Currently preemptible instances only supported for Google LifeSciences Executor
    if preemption_default or preemptible_rules and not google_lifesciences:
        logger.warning(
            "Preemptible instances are only available for the Google Life Sciences Executor."
        )

    if updated_files is None:
        updated_files = list()

    run_local = not (
        cluster
        or cluster_sync
        or drmaa
        or kubernetes
        or tibanna
        or az_batch
        or google_lifesciences
        or tes
        or slurm
        or slurm_jobstep
    )
    if run_local:
        if not dryrun:
            # clean up all previously recorded jobids.
            shell.cleanup()
    else:
        if default_resources is None:
            # use full default resources if in cluster or cloud mode
            default_resources = DefaultResources(mode="full")
        if edit_notebook:
            raise WorkflowError(
                "Notebook edit mode is only allowed with local execution."
            )

    shell.conda_block_conflicting_envvars = not conda_not_block_search_path_envvars

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
            printshellcmds=printshellcmds,
            debug_dag=debug_dag,
            nocolor=nocolor,
            stdout=stdout,
            debug=verbose,
            use_threads=use_threads,
            mode=mode,
            show_failed_logs=show_failed_logs,
            dryrun=dryrun,
        )

    if greediness is None:
        greediness = 0.5 if prioritytargets else 1.0
    else:
        if not (0 <= greediness <= 1.0):
            logger.error("Error: greediness must be a float between 0 and 1.")
            return False

    if not os.path.exists(snakefile):
        logger.error(f'Error: Snakefile "{snakefile}" not found.')
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
        update_config(overwrite_config, load_configfile(f))
    # convert provided paths to absolute paths
    configfiles = list(map(os.path.abspath, configfiles))

    # directly specified elements override any configfiles
    if config:
        update_config(overwrite_config, config)
        if config_args is None:
            config_args = dict_to_key_value_args(config)

    if workdir:
        olddir = os.getcwd()
        if not os.path.exists(workdir):
            logger.info(f"Creating specified working directory {workdir}.")
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
                    keep_local=keep_remote_local, is_default=True
                )
            else:
                raise WorkflowError(
                    "Remote provider {} does not (yet) support to "
                    "be used as default provider."
                )

        workflow = Workflow(
            snakefile=snakefile,
            rerun_triggers=rerun_triggers,
            jobscript=jobscript,
            overwrite_shellcmd=overwrite_shellcmd,
            overwrite_config=overwrite_config,
            overwrite_workdir=workdir,
            overwrite_configfiles=configfiles,
            overwrite_threads=overwrite_threads,
            max_threads=max_threads,
            overwrite_scatter=overwrite_scatter,
            overwrite_groups=overwrite_groups,
            overwrite_resources=overwrite_resources,
            overwrite_resource_scopes=overwrite_resource_scopes,
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
            assume_shared_fs=assume_shared_fs,
            default_resources=default_resources,
            cache=cache,
            cores=cores,
            nodes=nodes,
            resources=resources,
            edit_notebook=edit_notebook,
            envvars=envvars,
            max_inventory_wait_time=max_inventory_wait_time,
            conda_not_block_search_path_envvars=conda_not_block_search_path_envvars,
            execute_subworkflows=execute_subworkflows,
            scheduler_solver_path=scheduler_solver_path,
            conda_base_path=conda_base_path,
            check_envvars=not lint,  # for linting, we do not need to check whether requested envvars exist
            all_temp=all_temp,
            local_groupid=local_groupid,
            keep_metadata=keep_metadata,
            latency_wait=latency_wait,
            executor_args=executor_args,
            cleanup_scripts=cleanup_scripts,
            immediate_submit=immediate_submit,
            quiet=quiet,
        )
        success = True

        workflow.include(
            snakefile,
            overwrite_default_target=True,
            print_compilation=print_compilation,
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
                    max_threads=max_threads,
                    cache=cache,
                    overwrite_threads=overwrite_threads,
                    overwrite_scatter=overwrite_scatter,
                    overwrite_resources=overwrite_resources,
                    overwrite_resource_scopes=overwrite_resource_scopes,
                    default_resources=default_resources,
                    dryrun=dryrun,
                    touch=touch,
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
                    cleanup_containers=cleanup_containers,
                    cleanup_shadow=cleanup_shadow,
                    cleanup_scripts=cleanup_scripts,
                    force_incomplete=force_incomplete,
                    ignore_incomplete=ignore_incomplete,
                    latency_wait=latency_wait,
                    verbose=verbose,
                    notemp=notemp,
                    all_temp=all_temp,
                    keep_remote_local=keep_remote_local,
                    nodeps=nodeps,
                    jobscript=jobscript,
                    greediness=greediness,
                    no_hooks=no_hooks,
                    overwrite_shellcmd=overwrite_shellcmd,
                    config=config,
                    config_args=config_args,
                    keep_logger=True,
                    force_use_threads=use_threads,
                    use_conda=use_conda,
                    use_singularity=use_singularity,
                    use_env_modules=use_env_modules,
                    conda_prefix=conda_prefix,
                    conda_cleanup_pkgs=conda_cleanup_pkgs,
                    conda_frontend=conda_frontend,
                    singularity_prefix=singularity_prefix,
                    shadow_prefix=shadow_prefix,
                    singularity_args=singularity_args,
                    scheduler=scheduler,
                    scheduler_ilp_solver=scheduler_ilp_solver,
                    list_conda_envs=list_conda_envs,
                    kubernetes=kubernetes,
                    container_image=container_image,
                    k8s_cpu_scalar=k8s_cpu_scalar,
                    k8s_service_account_name=k8s_service_account_name,
                    conda_create_envs_only=conda_create_envs_only,
                    default_remote_provider=default_remote_provider,
                    default_remote_prefix=default_remote_prefix,
                    tibanna=tibanna,
                    tibanna_sfn=tibanna_sfn,
                    az_batch=az_batch,
                    az_batch_enable_autoscale=az_batch_enable_autoscale,
                    az_batch_account_url=az_batch_account_url,
                    google_lifesciences=google_lifesciences,
                    google_lifesciences_regions=google_lifesciences_regions,
                    google_lifesciences_location=google_lifesciences_location,
                    google_lifesciences_cache=google_lifesciences_cache,
                    google_lifesciences_service_account_email=google_lifesciences_service_account_email,
                    google_lifesciences_network=google_lifesciences_network,
                    google_lifesciences_subnetwork=google_lifesciences_subnetwork,
                    flux=flux,
                    tes=tes,
                    precommand=precommand,
                    preemption_default=preemption_default,
                    preemptible_rules=preemptible_rules,
                    tibanna_config=tibanna_config,
                    assume_shared_fs=assume_shared_fs,
                    cluster_status=cluster_status,
                    cluster_cancel=cluster_cancel,
                    cluster_cancel_nargs=cluster_cancel_nargs,
                    cluster_sidecar=cluster_sidecar,
                    max_jobs_per_second=max_jobs_per_second,
                    max_status_checks_per_second=max_status_checks_per_second,
                    overwrite_groups=overwrite_groups,
                    group_components=group_components,
                    max_inventory_wait_time=max_inventory_wait_time,
                    conda_not_block_search_path_envvars=conda_not_block_search_path_envvars,
                    local_groupid=local_groupid,
                )
                success = workflow.execute(
                    targets=targets,
                    target_jobs=target_jobs,
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
                    keepgoing=keepgoing,
                    printrulegraph=printrulegraph,
                    printfilegraph=printfilegraph,
                    printdag=printdag,
                    slurm=slurm,
                    slurm_jobstep=slurm_jobstep,
                    cluster=cluster,
                    cluster_sync=cluster_sync,
                    jobname=jobname,
                    drmaa=drmaa,
                    drmaa_log_dir=drmaa_log_dir,
                    kubernetes=kubernetes,
                    container_image=container_image,
                    k8s_cpu_scalar=k8s_cpu_scalar,
                    k8s_service_account_name=k8s_service_account_name,
                    tibanna=tibanna,
                    tibanna_sfn=tibanna_sfn,
                    az_batch=az_batch,
                    az_batch_enable_autoscale=az_batch_enable_autoscale,
                    az_batch_account_url=az_batch_account_url,
                    google_lifesciences=google_lifesciences,
                    google_lifesciences_regions=google_lifesciences_regions,
                    google_lifesciences_location=google_lifesciences_location,
                    google_lifesciences_cache=google_lifesciences_cache,
                    google_lifesciences_service_account_email=google_lifesciences_service_account_email,
                    google_lifesciences_network=google_lifesciences_network,
                    google_lifesciences_subnetwork=google_lifesciences_subnetwork,
                    tes=tes,
                    flux=flux,
                    precommand=precommand,
                    preemption_default=preemption_default,
                    preemptible_rules=preemptible_rules,
                    tibanna_config=tibanna_config,
                    max_jobs_per_second=max_jobs_per_second,
                    max_status_checks_per_second=max_status_checks_per_second,
                    printd3dag=printd3dag,
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
                    cleanup_containers=cleanup_containers,
                    cleanup_shadow=cleanup_shadow,
                    subsnakemake=subsnakemake,
                    updated_files=updated_files,
                    allowed_rules=allowed_rules,
                    greediness=greediness,
                    no_hooks=no_hooks,
                    force_use_threads=use_threads,
                    conda_create_envs_only=conda_create_envs_only,
                    cluster_status=cluster_status,
                    cluster_cancel=cluster_cancel,
                    cluster_cancel_nargs=cluster_cancel_nargs,
                    cluster_sidecar=cluster_sidecar,
                    report=report,
                    report_stylesheet=report_stylesheet,
                    export_cwl=export_cwl,
                    batch=batch,
                    keepincomplete=keep_incomplete,
                    containerize=containerize,
                )

    except BrokenPipeError:
        # ignore this exception and stop. It occurs if snakemake output is piped into less and less quits before reading the whole output.
        # in such a case, snakemake shall stop scheduling and quit with error 1
        success = False
    except BaseException as ex:
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
