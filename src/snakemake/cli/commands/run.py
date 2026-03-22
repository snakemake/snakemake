"""Click subcommand for executing a Snakemake workflow."""

from pathlib import Path
from typing import Optional, Tuple

import click

from snakemake.api import SnakemakeApi
from snakemake.cli.click import ProfileAwareCommand
from snakemake.cli.common import (
    workflow_options,
    resource_options,
    config_options,
    execution_options,
    dag_options,
    output_options,
    deployment_options,
    storage_options,
    scheduling_options,
    remote_execution_options,
    group_options,
)
from snakemake.cli.legacy import (
    parse_batch,
    parse_config,
    parse_consider_ancient,
    parse_groups,
    parse_group_components,
    parse_set_resources,
    parse_set_scatter,
    parse_set_threads,
    parse_set_resource_scope,
    parse_size_in_bytes,
    parse_timespan,
)
from snakemake.cli.operations import run_report
from snakemake.cli.plugin_adapter import add_plugin_options, plugin_settings_from_kwargs
from snakemake.common import get_container_image
from snakemake.settings.enums import (
    CondaCleanupPkgs,
    Quietness,
    RerunTrigger,
    StrictDagEvaluation,
)
from snakemake.settings.types import (
    ConfigSettings,
    DAGSettings,
    DeploymentSettings,
    ExecutionSettings,
    GroupSettings,
    MaxJobsPerTimespan,
    OutputSettings,
    PreemptibleRules,
    RemoteExecutionSettings,
    ResourceSettings,
    SchedulingSettings,
    StorageSettings,
    WorkflowSettings,
)
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
from snakemake_interface_executor_plugins.settings import (
    DeploymentMethod,
    SharedFSUsage,
)
from snakemake_interface_logger_plugins.registry import LoggerPluginRegistry
from snakemake_interface_report_plugins.registry import ReportPluginRegistry
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry
from snakemake_interface_scheduler_plugins.registry import SchedulerPluginRegistry
from snakemake.resources import Resources


@click.command(cls=ProfileAwareCommand)
@workflow_options
@resource_options
@config_options
@execution_options
@dag_options
@output_options
@deployment_options
@storage_options
@scheduling_options
@remote_execution_options
@group_options
@click.option(
    "-e",
    "--executor",
    type=click.Choice(list(ExecutorPluginRegistry().plugins.keys())),
    default=None,
    help="Specify a custom executor plugin.",
)
@click.option(
    "--scheduler",
    type=click.Choice(list(SchedulerPluginRegistry().plugins.keys())),
    default="ilp",
    show_default=True,
    help="Scheduling algorithm.",
)
@click.option(
    "--logger",
    "loggers",
    multiple=True,
    type=click.Choice(list(LoggerPluginRegistry().plugins.keys())),
    help="Specify one or more logger plugins.",
)
@click.option(
    "--reporter",
    type=click.Choice(list(ReportPluginRegistry().plugins.keys())),
    default=None,
    help="Report plugin to use.",
)
@click.option(
    "--profile",
    default=None,
    envvar="SNAKEMAKE_PROFILE",
    help="Profile to use for setting default options.",
)
@click.option(
    "--workflow-profile",
    default=None,
    help="Workflow-specific profile for setting default options.",
)
@click.option(
    "--report-after-run",
    is_flag=True,
    default=False,
    help="Generate a report after workflow execution.",
)
@click.option(
    "--report",
    type=click.Path(path_type=Path),
    default=None,
    help="Report output path (required with --report-after-run).",
)
@click.option(
    "--report-stylesheet",
    type=click.Path(path_type=Path),
    default=None,
    help="Custom stylesheet for the report.",
)
@click.option(
    "--report-metadata",
    type=click.Path(path_type=Path),
    default=None,
    help="Custom metadata file for the report.",
)
@click.option(
    "--cache",
    multiple=True,
    metavar="RULE",
    help="Cache output files of given rules.",
)
@click.option(
    "--consider-ancient",
    multiple=True,
    metavar="RULE=INPUTITEMS",
    help="Consider given input items of given rules as ancient.",
)
@click.option(
    "--wrapper-prefix",
    default="https://github.com/snakemake/snakemake-wrappers/raw/",
    show_default=True,
    help="URL prefix for the wrapper directive.",
)
@add_plugin_options("executor", "storage", "report", "logger", "scheduler")
def run(
    # workflow_options
    snakefile: Optional[Path],
    directory: Optional[Path],
    # resource_options
    cores: Optional[int],
    jobs: Optional[int],
    local_cores: Optional[int],
    max_threads: Optional[int],
    set_threads: Tuple[str, ...],
    set_resources: Tuple[str, ...],
    set_scatter: Tuple[str, ...],
    set_resource_scopes: Tuple[str, ...],
    default_resources: Tuple[str, ...],
    resources: Tuple[str, ...],
    # config_options
    configfile: Tuple[Path, ...],
    config: Tuple[str, ...],
    replace_workflow_config: bool,
    # execution_options
    dryrun: bool,
    touch: bool,
    keep_going: bool,
    debug: bool,
    latency_wait: int,
    retries: int,
    nolock: bool,
    ignore_incomplete: bool,
    rerun_incomplete: bool,
    shadow_prefix: Optional[Path],
    force_use_threads: bool,
    allow_ambiguity: bool,
    keep_incomplete: bool,
    drop_metadata: bool,
    skip_script_cleanup: bool,
    no_hooks: bool,
    queue_input_wait_time: int,
    show_failed_logs: bool,
    # dag_options
    targets: Tuple[str, ...],
    force: bool,
    forceall: bool,
    forcerun: Tuple[str, ...],
    until: Tuple[str, ...],
    omit_from: Tuple[str, ...],
    batch: Optional[str],
    rerun_triggers: Tuple[str, ...],
    allowed_rules: Tuple[str, ...],
    max_inventory_time: int,
    trust_io_cache: bool,
    max_checksum_file_size: str,
    strict_dag_evaluation: Tuple[str, ...],
    # output_options
    printshellcmds: bool,
    nocolor: bool,
    quiet: Tuple[str, ...],
    verbose: bool,
    debug_dag: bool,
    benchmark_extended: bool,
    # deployment_options
    software_deployment_method: Tuple[str, ...],
    conda_prefix: Optional[Path],
    conda_frontend: str,
    conda_not_block_search_path_envvars: bool,
    conda_base_path: Optional[Path],
    conda_cleanup_pkgs: Optional[str],
    apptainer_args: str,
    apptainer_prefix: Optional[Path],
    # storage_options
    default_storage_provider: Optional[str],
    default_storage_prefix: str,
    local_storage_prefix: Path,
    remote_job_local_storage_prefix: Path,
    shared_fs_usage: Tuple[str, ...],
    notemp: bool,
    all_temp: bool,
    keep_storage_local_copies: bool,
    not_retrieve_storage: bool,
    omit_flags: Tuple[str, ...],
    unneeded_temp_files: Tuple[str, ...],
    wait_for_free_local_storage: Optional[str],
    # scheduling_options (--scheduler is below, not in the decorator)
    prioritize: Tuple[str, ...],
    scheduler_greediness: Optional[float],
    scheduler_subsample: Optional[int],
    max_jobs_per_timespan: str,
    # remote_execution_options
    container_image: Optional[str],
    jobname: str,
    jobscript: Optional[Path],
    envvars: Tuple[str, ...],
    immediate_submit: bool,
    max_status_checks_per_second: float,
    seconds_between_status_checks: int,
    precommand: Optional[str],
    job_deploy_sources: bool,
    preemptible_rules: Tuple[str, ...],
    preemptible_retries: Optional[int],
    # group_options
    groups: Tuple[str, ...],
    group_components: Tuple[str, ...],
    # run-specific options
    executor: Optional[str],
    scheduler: str,
    loggers: Tuple[str, ...],
    reporter: Optional[str],
    profile: Optional[str],
    workflow_profile: Optional[str],
    report_after_run: bool,
    report: Optional[Path],
    report_stylesheet: Optional[Path],
    report_metadata: Optional[Path],
    cache: Tuple[str, ...],
    consider_ancient: Tuple[str, ...],
    wrapper_prefix: str,
    # plugin options (dynamically added)
    **plugin_kwargs,
):
    """Execute the workflow."""

    executor_name = executor
    if dryrun:
        executor_name = "dryrun"
    elif touch:
        executor_name = "touch"
    elif executor_name is None:
        executor_name = "local"

    executor_plugin = ExecutorPluginRegistry().get_plugin(executor_name)
    executor_settings = plugin_settings_from_kwargs(executor_plugin, plugin_kwargs)

    if cores is None:
        if executor_plugin.common_settings.local_exec:  # type: ignore[attr-defined]
            cores = jobs
            jobs = None
        elif executor_plugin.common_settings.dryrun_exec:  # type: ignore[attr-defined]
            cores = 1
            jobs = None

    storage_provider_settings = {
        name: plugin_settings_from_kwargs(
            StoragePluginRegistry().get_plugin(name), plugin_kwargs
        )
        for name in StoragePluginRegistry().get_registered_plugins()
    }

    log_handler_settings = {
        name: plugin_settings_from_kwargs(
            LoggerPluginRegistry().get_plugin(name), plugin_kwargs
        )
        for name in loggers
    }

    scheduler_plugin = SchedulerPluginRegistry().get_plugin(scheduler)
    scheduler_settings = plugin_settings_from_kwargs(scheduler_plugin, plugin_kwargs)

    if reporter:
        report_plugin = ReportPluginRegistry().get_plugin(reporter)
        report_plugin_settings = plugin_settings_from_kwargs(
            report_plugin, plugin_kwargs
        )
    elif report_after_run:
        raise click.UsageError("--report-after-run requires --reporter to be set.")
    else:
        report_plugin = None
        report_plugin_settings = None

    parsed_rerun_triggers = RerunTrigger.parse_choices_set(list(rerun_triggers))
    parsed_deployment_method = DeploymentMethod.parse_choices_set(
        list(software_deployment_method)
    )
    parsed_shared_fs_usage = SharedFSUsage.parse_choices_set(list(shared_fs_usage))
    parsed_quiet = Quietness.parse_choices_set(list(quiet)) if quiet else None
    parsed_max_jobs_per_timespan = MaxJobsPerTimespan.parse_choice(
        max_jobs_per_timespan
    )
    parsed_batch = parse_batch(batch) if batch else None
    parsed_config = parse_config(config) if config else {}
    parsed_set_threads = parse_set_threads(list(set_threads)) if set_threads else {}
    parsed_set_resources = (
        parse_set_resources(list(set_resources)) if set_resources else {}
    )
    parsed_set_scatter = parse_set_scatter(list(set_scatter)) if set_scatter else {}
    parsed_set_resource_scopes = (
        parse_set_resource_scope(list(set_resource_scopes))
        if set_resource_scopes
        else {}
    )
    parsed_groups = parse_groups(list(groups)) if groups else {}
    parsed_group_components = (
        parse_group_components(list(group_components)) if group_components else {}
    )
    parsed_consider_ancient = (
        parse_consider_ancient(list(consider_ancient)) if consider_ancient else {}
    )
    parsed_strict_dag = StrictDagEvaluation.parse_choices_set(
        list(strict_dag_evaluation)
    )
    parsed_checksum_size = parse_size_in_bytes(max_checksum_file_size)
    parsed_conda_cleanup = (
        CondaCleanupPkgs.parse_choice(conda_cleanup_pkgs)
        if conda_cleanup_pkgs
        else None
    )
    parsed_wait_storage = (
        parse_timespan(wait_for_free_local_storage)
        if wait_for_free_local_storage
        else None
    )

    parsed_resources = (
        Resources.parse(list(resources), allow_expressions=False)
        if resources
        else Resources()
    )
    parsed_default_resources = (
        Resources.parse(
            list(default_resources), defaults="full", allow_expressions=True
        )
        if default_resources
        else None
    )
    if len(preemptible_rules) > 0:
        parsed_preemptible = PreemptibleRules(rules=set(preemptible_rules))
    else:
        parsed_preemptible = PreemptibleRules()

    output_settings = OutputSettings(
        dryrun=dryrun,
        printshellcmds=printshellcmds,
        nocolor=nocolor,
        quiet=parsed_quiet,
        debug_dag=debug_dag,
        verbose=verbose,
        show_failed_logs=show_failed_logs,
        log_handler_settings=log_handler_settings,
        stdout=dryrun,
        benchmark_extended=benchmark_extended,
    )

    resource_settings = ResourceSettings(
        cores=cores,
        nodes=jobs,
        local_cores=local_cores,
        max_threads=max_threads,
        resources=parsed_resources,
        default_resources=parsed_default_resources,
        overwrite_threads=parsed_set_threads,
        overwrite_scatter=parsed_set_scatter,
        overwrite_resource_scopes=parsed_set_resource_scopes,
        overwrite_resources=parsed_set_resources,
    )

    config_settings = ConfigSettings(
        config=parsed_config,
        configfiles=list(configfile) if configfile else [],
        config_args=None,
        replace_workflow_config=replace_workflow_config,
    )

    storage_settings = StorageSettings(
        default_storage_provider=default_storage_provider,
        default_storage_prefix=default_storage_prefix,
        local_storage_prefix=local_storage_prefix,
        remote_job_local_storage_prefix=remote_job_local_storage_prefix,
        shared_fs_usage=parsed_shared_fs_usage,
        keep_storage_local=keep_storage_local_copies,
        retrieve_storage=not not_retrieve_storage,
        notemp=notemp,
        all_temp=all_temp,
        omit_flags=set(omit_flags) if omit_flags else frozenset(),
        unneeded_temp_files=set(unneeded_temp_files)
        if unneeded_temp_files
        else frozenset(),
        wait_for_free_local_storage=parsed_wait_storage,  # type: ignore
    )

    workflow_settings = WorkflowSettings(
        wrapper_prefix=wrapper_prefix,
        cache=list(cache) if cache else None,
        consider_ancient=parsed_consider_ancient,
    )

    deployment_settings = DeploymentSettings(
        deployment_method=parsed_deployment_method,
        conda_prefix=conda_prefix,
        conda_frontend=conda_frontend,
        conda_not_block_search_path_envvars=conda_not_block_search_path_envvars,
        conda_base_path=conda_base_path,
        conda_cleanup_pkgs=parsed_conda_cleanup,  # type: ignore
        apptainer_args=apptainer_args,
        apptainer_prefix=apptainer_prefix,
    )

    dag_settings = DAGSettings(
        targets=set(targets) if targets else set(),
        forcetargets=force,
        forceall=forceall,
        forcerun=set(forcerun) if forcerun else set(),
        until=set(until) if until else set(),
        omit_from=set(omit_from) if omit_from else set(),
        batch=parsed_batch,
        force_incomplete=rerun_incomplete,
        rerun_triggers=parsed_rerun_triggers,
        allowed_rules=set(allowed_rules) if allowed_rules else set(),
        max_inventory_wait_time=max_inventory_time,
        trust_io_cache=trust_io_cache,
        max_checksum_file_size=parsed_checksum_size,
        strict_evaluation=parsed_strict_dag,
    )

    execution_settings = ExecutionSettings(
        keep_going=keep_going,
        debug=debug,
        standalone=True,
        ignore_ambiguity=allow_ambiguity,
        lock=not nolock,
        ignore_incomplete=ignore_incomplete,
        latency_wait=latency_wait,
        no_hooks=no_hooks,
        retries=retries,
        use_threads=force_use_threads,
        shadow_prefix=shadow_prefix,
        keep_incomplete=keep_incomplete,
        keep_metadata=not drop_metadata,
        cleanup_scripts=not skip_script_cleanup,
        queue_input_wait_time=queue_input_wait_time,
    )

    remote_execution_settings = RemoteExecutionSettings(
        jobname=jobname,
        jobscript=jobscript,
        max_status_checks_per_second=max_status_checks_per_second,
        seconds_between_status_checks=seconds_between_status_checks,
        container_image=container_image or get_container_image(),
        envvars=envvars,
        immediate_submit=immediate_submit,
        precommand=precommand,
        job_deploy_sources=job_deploy_sources,
        preemptible_retries=preemptible_retries,
        preemptible_rules=parsed_preemptible,
    )

    scheduling_settings = SchedulingSettings(
        prioritytargets=set(prioritize) if prioritize else set(),
        scheduler=scheduler,
        greediness=scheduler_greediness,
        subsample=scheduler_subsample,
        max_jobs_per_timespan=parsed_max_jobs_per_timespan,
    )

    group_settings = GroupSettings(
        overwrite_groups=parsed_groups,
        group_components=parsed_group_components,
    )

    with SnakemakeApi(output_settings) as snakemake_api:
        try:
            workflow_api = snakemake_api.workflow(
                resource_settings=resource_settings,
                config_settings=config_settings,
                storage_settings=storage_settings,
                storage_provider_settings=storage_provider_settings,
                workflow_settings=workflow_settings,
                deployment_settings=deployment_settings,
                snakefile=snakefile,
                workdir=directory,
            )

            dag_api = workflow_api.dag(dag_settings=dag_settings)

            dag_api.execute_workflow(
                executor=executor_name,
                execution_settings=execution_settings,
                remote_execution_settings=remote_execution_settings,
                scheduling_settings=scheduling_settings,
                group_settings=group_settings,
                executor_settings=executor_settings,
                scheduler_settings=scheduler_settings,
            )

            if report_plugin is not None and report_after_run:
                run_report(
                    dag_api,
                    reporter,  # type: ignore
                    report_plugin_settings,  # type: ignore
                    report_metadata,
                )

        except Exception as e:
            snakemake_api.print_exception(e)
            raise SystemExit(1)
