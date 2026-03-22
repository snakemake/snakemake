"""Shared options and context for Click subcommands.

Each decorator group corresponds to a settings dataclass in
``snakemake.settings.types``. The ``run`` command stacks all of them;
simpler subcommands use only ``workflow_options``.
"""

import functools
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Optional

import click

from snakemake.api import SnakemakeApi
from snakemake.settings.types import (
    ConfigSettings,
    DAGSettings,
    DeploymentSettings,
    OutputSettings,
    ResourceSettings,
    StorageSettings,
    WorkflowSettings,
)


def workflow_options(func):
    """Core options shared by every subcommand: --snakefile, --directory."""

    @click.option(
        "-s",
        "--snakefile",
        type=click.Path(path_type=Path),
        default=None,
        help="Path to the Snakefile.",
    )
    @click.option(
        "-d",
        "--directory",
        type=click.Path(path_type=Path),
        default=None,
        help="Specify working directory.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def resource_options(func):
    """Options that feed into ``ResourceSettings``."""

    @click.option(
        "-c",
        "--cores",
        type=int,
        default=None,
        help="Use at most N CPU cores/jobs in parallel.",
    )
    @click.option(
        "-j",
        "--jobs",
        type=int,
        default=None,
        help="Use at most N CPU cluster/cloud jobs in parallel.",
    )
    @click.option(
        "--local-cores",
        type=int,
        default=None,
        help="In cluster/cloud mode, use at most N cores for local rules.",
    )
    @click.option(
        "--max-threads",
        type=int,
        default=None,
        help="Global maximum number of threads available to any rule.",
    )
    @click.option(
        "--set-threads",
        multiple=True,
        metavar="RULE=THREADS",
        help="Overwrite thread usage of rules.",
    )
    @click.option(
        "--set-resources",
        multiple=True,
        metavar="RULE:RESOURCE=VALUE",
        help="Overwrite resource usage of rules.",
    )
    @click.option(
        "--set-scatter",
        multiple=True,
        metavar="NAME=SCATTERITEMS",
        help="Overwrite number of scatter items.",
    )
    @click.option(
        "--set-resource-scopes",
        multiple=True,
        metavar="RESOURCE=[global|local]",
        help="Overwrite resource scopes.",
    )
    @click.option(
        "--default-resources",
        multiple=True,
        metavar="NAME=INT",
        help="Define default values of resources for rules.",
    )
    @click.option(
        "--resources",
        "--res",
        multiple=True,
        metavar="NAME=INT",
        help="Define additional resources that constrain scheduling (e.g. mem_mb=1000).",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def config_options(func):
    """Options that feed into ``ConfigSettings``."""

    @click.option(
        "--configfile",
        "--configfiles",
        multiple=True,
        type=click.Path(path_type=Path),
        help="Specify or overwrite the config file of the workflow.",
    )
    @click.option(
        "-C",
        "--config",
        multiple=True,
        metavar="KEY=VALUE",
        help="Set or overwrite values in the workflow config object.",
    )
    @click.option(
        "--replace-workflow-config",
        is_flag=True,
        default=False,
        help="Config files fully replace the workflow config instead of updating it.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def execution_options(func):
    """Options that feed into ``ExecutionSettings``."""

    @click.option(
        "-n",
        "--dry-run",
        "--dryrun",
        "dryrun",
        is_flag=True,
        default=False,
        help="Do not execute anything, only display what would be done.",
    )
    @click.option(
        "-t",
        "--touch",
        is_flag=True,
        default=False,
        help="Touch output files (mark them up to date) instead of running.",
    )
    @click.option(
        "-k",
        "--keep-going",
        is_flag=True,
        default=False,
        help="Go on with independent jobs if a job fails.",
    )
    @click.option(
        "--debug",
        is_flag=True,
        default=False,
        help="Allow debugging rules with e.g. PDB.",
    )
    @click.option(
        "-w",
        "--latency-wait",
        type=int,
        default=5,
        show_default=True,
        help="Wait given seconds if an output file is not present after the job finished.",
    )
    @click.option(
        "-T",
        "--retries",
        type=int,
        default=0,
        show_default=True,
        help="Number of times to restart failing jobs.",
    )
    @click.option(
        "--nolock",
        is_flag=True,
        default=False,
        help="Do not lock the working directory.",
    )
    @click.option(
        "--ignore-incomplete",
        "--ii",
        is_flag=True,
        default=False,
        help="Do not check for incomplete output files.",
    )
    @click.option(
        "--rerun-incomplete",
        "--ri",
        is_flag=True,
        default=False,
        help="Re-run all jobs the output of which is recognized as incomplete.",
    )
    @click.option(
        "--shadow-prefix",
        type=click.Path(path_type=Path),
        default=None,
        help="Specify a directory in which the shadow directory is created.",
    )
    @click.option(
        "--force-use-threads",
        is_flag=True,
        default=False,
        help="Force threads rather than processes.",
    )
    @click.option(
        "--allow-ambiguity",
        "-a",
        is_flag=True,
        default=False,
        help="Don't check for ambiguous rules.",
    )
    @click.option(
        "--keep-incomplete",
        is_flag=True,
        default=False,
        help="Do not remove incomplete output files by failed jobs.",
    )
    @click.option(
        "--drop-metadata",
        is_flag=True,
        default=False,
        help="Drop metadata file tracking information after job finishes.",
    )
    @click.option(
        "--skip-script-cleanup",
        is_flag=True,
        default=False,
        help="Don't delete wrapper scripts used for execution.",
    )
    @click.option(
        "--no-hooks",
        is_flag=True,
        default=False,
        help="Do not invoke onstart, onsuccess or onerror hooks.",
    )
    @click.option(
        "--queue-input-wait-time",
        type=int,
        default=10,
        show_default=True,
        help="Interval in seconds to check for new input in queue rules.",
    )
    @click.option(
        "--show-failed-logs",
        is_flag=True,
        default=False,
        help="Automatically display logs of failed jobs.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def dag_options(func):
    """Options that feed into ``DAGSettings``."""

    @click.argument("targets", nargs=-1)
    @click.option(
        "-f",
        "--force",
        is_flag=True,
        default=False,
        help="Force the execution of the selected target or the first rule.",
    )
    @click.option(
        "-F",
        "--forceall",
        is_flag=True,
        default=False,
        help="Force execution of all rules regardless of existing output.",
    )
    @click.option(
        "-R",
        "--forcerun",
        multiple=True,
        metavar="TARGET",
        help="Force the re-execution of the given rules or files.",
    )
    @click.option(
        "-U",
        "--until",
        multiple=True,
        metavar="TARGET",
        help="Run the pipeline until it reaches the specified rules or files.",
    )
    @click.option(
        "-O",
        "--omit-from",
        multiple=True,
        metavar="TARGET",
        help="Prevent execution of the given rules and their downstream dependencies.",
    )
    @click.option(
        "--batch",
        default=None,
        metavar="RULE=BATCH/BATCHES",
        help="Only create the given BATCH of input files of the given RULE.",
    )
    @click.option(
        "--rerun-triggers",
        multiple=True,
        type=click.Choice(["mtime", "params", "input", "software-env", "code"]),
        default=("mtime", "params", "input", "software-env", "code"),
        help="Define what triggers the rerunning of a job.",
    )
    @click.option(
        "--allowed-rules",
        multiple=True,
        help="Only consider given rules.",
    )
    @click.option(
        "--max-inventory-time",
        type=int,
        default=20,
        show_default=True,
        help="Spend at most SECONDS creating a file inventory.",
    )
    @click.option(
        "--trust-io-cache",
        is_flag=True,
        default=False,
        help="Assume all input/output queries from previous dryruns are still valid.",
    )
    @click.option(
        "--max-checksum-file-size",
        default="1MB",
        help="Only checksum files smaller than this threshold (e.g. 1MB).",
    )
    @click.option(
        "--strict-dag-evaluation",
        multiple=True,
        type=click.Choice(["functions", "cyclic-graph", "periodic-wildcards"]),
        default=(),
        help="Strict evaluation of rules even when not required for output.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def output_options(func):
    """Options that feed into ``OutputSettings``."""

    @click.option(
        "-p",
        "--printshellcmds",
        is_flag=True,
        default=False,
        help="Print out the shell commands that will be executed.",
    )
    @click.option(
        "--nocolor",
        is_flag=True,
        default=False,
        help="Do not use a colored output.",
    )
    @click.option(
        "-q",
        "--quiet",
        multiple=True,
        type=click.Choice(["rules", "progress", "all", "host", "reason"]),
        default=None,
        help="Do not output certain information.",
    )
    @click.option(
        "--verbose",
        is_flag=True,
        default=False,
        help="Print debugging output.",
    )
    @click.option(
        "--debug-dag",
        is_flag=True,
        default=False,
        help="Print candidate and selected jobs while inferring DAG.",
    )
    @click.option(
        "--benchmark-extended",
        is_flag=True,
        default=False,
        help="Write extended benchmarking metrics.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def deployment_options(func):
    """Options that feed into ``DeploymentSettings``."""

    @click.option(
        "--software-deployment-method",
        "--sdm",
        multiple=True,
        type=click.Choice(["conda", "apptainer", "env-modules"]),
        default=(),
        help="Specify software environment deployment method.",
    )
    @click.option(
        "--conda-prefix",
        type=click.Path(path_type=Path),
        default=os.environ.get("SNAKEMAKE_CONDA_PREFIX", None),
        envvar="SNAKEMAKE_CONDA_PREFIX",
        help="Directory in which conda environments are stored.",
    )
    @click.option(
        "--conda-frontend",
        type=click.Choice(["conda", "mamba"]),
        default="conda",
        show_default=True,
        help="Choose the conda frontend.",
    )
    @click.option(
        "--conda-not-block-search-path-envvars",
        is_flag=True,
        default=False,
        help="Do not block search path environment variables when using conda.",
    )
    @click.option(
        "--conda-base-path",
        type=click.Path(path_type=Path),
        default=None,
        help="Path to conda base installation.",
    )
    @click.option(
        "--conda-cleanup-pkgs",
        type=click.Choice(["tarballs", "cache"]),
        default=None,
        help="Cleanup conda packages after creating environments.",
    )
    @click.option(
        "--apptainer-args",
        default="",
        metavar="ARGS",
        help="Pass additional args to apptainer/singularity.",
    )
    @click.option(
        "--apptainer-prefix",
        type=click.Path(path_type=Path),
        default=None,
        help="Directory in which apptainer/singularity images are stored.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def storage_options(func):
    """Options that feed into ``StorageSettings``."""

    @click.option(
        "--default-storage-provider",
        default=None,
        help="Specify default storage provider for all input/output files.",
    )
    @click.option(
        "--default-storage-prefix",
        default="",
        help="Specify prefix for default storage provider.",
    )
    @click.option(
        "--local-storage-prefix",
        type=click.Path(path_type=Path),
        default=".snakemake/storage",
        show_default=True,
        help="Prefix for storing local copies of storage files.",
    )
    @click.option(
        "--remote-job-local-storage-prefix",
        type=click.Path(path_type=Path),
        default=".snakemake/storage",
        show_default=True,
        help="Prefix for local copies in remote jobs.",
    )
    @click.option(
        "--shared-fs-usage",
        multiple=True,
        type=click.Choice(
            [
                "input-output",
                "persistence",
                "software-deployment",
                "sources",
                "source-cache",
                "storage-local-copies",
                "none",
            ]
        ),
        default=None,
        help="Set assumptions on shared filesystem for non-local execution.",
    )
    @click.option(
        "--notemp",
        "--no-temp",
        "notemp",
        is_flag=True,
        default=False,
        help="Ignore temp() declarations.",
    )
    @click.option(
        "--all-temp",
        is_flag=True,
        default=False,
        help="Mark all output files as temp files.",
    )
    @click.option(
        "--keep-storage-local-copies",
        is_flag=True,
        default=False,
        help="Keep local copies of remote input and output files.",
    )
    @click.option(
        "--not-retrieve-storage",
        is_flag=True,
        default=False,
        help="Do not retrieve remote files.",
    )
    @click.option(
        "--omit-flags",
        multiple=True,
        help="Omit the given input and output file flags (e.g. pipe).",
    )
    @click.option(
        "--unneeded-temp-files",
        multiple=True,
        metavar="FILE",
        help="Files that will not be uploaded and are deleted immediately after job completion.",
    )
    @click.option(
        "--wait-for-free-local-storage",
        default=None,
        help="Wait for given timespan for enough free local storage (e.g. 1h, 30m).",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def scheduling_options(func):
    """Options that feed into ``SchedulingSettings``.

    Note: --scheduler is not included here because it needs a
    ``click.Choice`` populated from ``SchedulerPluginRegistry``.
    It lives directly on the ``run`` command instead.
    """

    @click.option(
        "-P",
        "--prioritize",
        multiple=True,
        metavar="TARGET",
        help="Assign highest priority to creation of given targets.",
    )
    @click.option(
        "--scheduler-greediness",
        type=float,
        default=None,
        help="Set the greediness of scheduling (0 to 1).",
    )
    @click.option(
        "--scheduler-subsample",
        type=int,
        default=None,
        help="Set the number of jobs to be considered for scheduling.",
    )
    @click.option(
        "--max-jobs-per-timespan",
        default="100/1s",
        help="Maximal number of submissions per timespan (e.g. 50/1m).",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def remote_execution_options(func):
    """Options that feed into ``RemoteExecutionSettings``."""

    @click.option(
        "--container-image",
        default=None,
        help="Docker image to use for remote execution.",
    )
    @click.option(
        "--jobname",
        default="snakejob.{name}.{jobid}.sh",
        show_default=True,
        help="Provide a custom name for the cluster jobscript.",
    )
    @click.option(
        "--jobscript",
        type=click.Path(path_type=Path),
        default=None,
        help="Provide a custom job script for cluster submission.",
    )
    @click.option(
        "--envvars",
        multiple=True,
        metavar="VARNAME",
        help="Environment variables to pass to cloud jobs.",
    )
    @click.option(
        "--immediate-submit",
        is_flag=True,
        default=False,
        help="Immediately submit all jobs to the cluster.",
    )
    @click.option(
        "--max-status-checks-per-second",
        type=float,
        default=100.0,
        show_default=True,
        help="Maximal number of job status checks per second.",
    )
    @click.option(
        "--seconds-between-status-checks",
        type=int,
        default=10,
        show_default=True,
        help="Seconds to wait between two rounds of status checks.",
    )
    @click.option(
        "--precommand",
        default=None,
        help="Command to execute before each remote job.",
    )
    @click.option(
        "--job-deploy-sources",
        is_flag=True,
        default=True,
        help="Deploy workflow sources before remote job starts.",
    )
    @click.option(
        "--preemptible-rules",
        multiple=True,
        metavar="RULE",
        help="Rules that may use preemptible machines. If given with no rules, all are preemptible.",
    )
    @click.option(
        "--preemptible-retries",
        type=int,
        default=None,
        help="Number of retries for preemptible job failures.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def group_options(func):
    """Options that feed into ``GroupSettings``."""

    @click.option(
        "--groups",
        multiple=True,
        metavar="RULE=GROUP",
        help="Assign rules to groups.",
    )
    @click.option(
        "--group-components",
        multiple=True,
        metavar="GROUP=COMPONENTS",
        help="Set the number of connected components a group may span.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


class SnakemakeContext:
    """Bridges click subcommands and the Snakemake API.

    Holds all settings needed to construct API objects. Context managers
    handle SnakemakeApi enter/exit so callers get proper cleanup.
    """

    def __init__(
        self,
        snakefile: Optional[Path] = None,
        workdir: Optional[Path] = None,
        output_settings: Optional[OutputSettings] = None,
        resource_settings: Optional[ResourceSettings] = None,
        config_settings: Optional[ConfigSettings] = None,
        storage_settings: Optional[StorageSettings] = None,
        workflow_settings: Optional[WorkflowSettings] = None,
        deployment_settings: Optional[DeploymentSettings] = None,
        dag_settings: Optional[DAGSettings] = None,
    ):
        self.snakefile = Path(snakefile) if snakefile else None
        self.workdir = Path(workdir) if workdir else None
        self.output_settings = output_settings or OutputSettings()
        self.resource_settings = resource_settings or ResourceSettings()
        self.config_settings = config_settings or ConfigSettings()
        self.storage_settings = storage_settings or StorageSettings()
        self.workflow_settings = workflow_settings or WorkflowSettings()
        self.deployment_settings = deployment_settings or DeploymentSettings()
        self.dag_settings = dag_settings or DAGSettings()

    @contextmanager
    def workflow(self):
        """Yield a WorkflowApi inside a managed SnakemakeApi context."""
        with SnakemakeApi(self.output_settings) as api:
            try:
                workflow_api = api.workflow(
                    resource_settings=self.resource_settings,
                    config_settings=self.config_settings,
                    storage_settings=self.storage_settings,
                    workflow_settings=self.workflow_settings,
                    deployment_settings=self.deployment_settings,
                    snakefile=self.snakefile,
                    workdir=self.workdir,
                )
                yield workflow_api
            except Exception as e:
                api.print_exception(e)
                raise SystemExit(1)

    @contextmanager
    def dag(self):
        """Yield a DAGApi inside a managed SnakemakeApi context."""
        with self.workflow() as workflow_api:
            yield workflow_api.dag(dag_settings=self.dag_settings)