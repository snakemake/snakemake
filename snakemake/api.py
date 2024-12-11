__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import ABC
from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import sys
from typing import Dict, List, Mapping, Optional, Set
import os
import tarfile

from snakemake.common import MIN_PY_VERSION, SNAKEFILE_CHOICES, async_run
from snakemake.settings.types import (
    ChangeType,
    GroupSettings,
    SchedulingSettings,
    WorkflowSettings,
)

if sys.version_info < MIN_PY_VERSION:
    raise ValueError(f"Snakemake requires at least Python {'.'.join(MIN_PY_VERSION)}.")

from snakemake.common.workdir_handler import WorkdirHandler
from snakemake.settings.types import (
    DAGSettings,
    DeploymentMethod,
    DeploymentSettings,
    ExecutionSettings,
    OutputSettings,
    ConfigSettings,
    RemoteExecutionSettings,
    ResourceSettings,
    StorageSettings,
    SharedFSUsage,
)

from snakemake_interface_executor_plugins.settings import ExecMode, ExecutorSettingsBase
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
from snakemake_interface_common.exceptions import ApiError
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry
from snakemake_interface_common.plugin_registry.plugin import TaggedSettings
from snakemake_interface_report_plugins.settings import ReportSettingsBase
from snakemake_interface_report_plugins.registry import ReportPluginRegistry

from snakemake.workflow import Workflow
from snakemake.exceptions import print_exception
from snakemake.logging import setup_logger, logger
from snakemake.shell import shell
from snakemake.common import (
    MIN_PY_VERSION,
    __version__,
)
from snakemake.resources import DefaultResources


class ApiBase(ABC):
    def __post_init__(self):
        self._check()

    def _check(self):
        # nothing to check by default
        # override in subclasses if needed
        pass


def resolve_snakefile(path: Optional[Path], allow_missing: bool = False):
    """Get path to the snakefile.

    Arguments
    ---------
    path: Optional[Path] -- The path to the snakefile. If not provided, default locations will be tried.
    """
    if path is None:
        for p in SNAKEFILE_CHOICES:
            if p.exists():
                return p
        if not allow_missing:
            raise ApiError(
                f"No Snakefile found, tried {', '.join(map(str, SNAKEFILE_CHOICES))}."
            )
    return path


@dataclass
class SnakemakeApi(ApiBase):
    """The Snakemake API.

    Arguments
    ---------

    output_settings: OutputSettings -- The output settings for the Snakemake API.
    """

    output_settings: OutputSettings = field(default_factory=OutputSettings)
    _workflow_api: Optional["WorkflowApi"] = field(init=False, default=None)
    _is_in_context: bool = field(init=False, default=False)

    def workflow(
        self,
        resource_settings: ResourceSettings,
        config_settings: Optional[ConfigSettings] = None,
        storage_settings: Optional[StorageSettings] = None,
        workflow_settings: Optional[WorkflowSettings] = None,
        deployment_settings: Optional[DeploymentSettings] = None,
        storage_provider_settings: Optional[Mapping[str, TaggedSettings]] = None,
        snakefile: Optional[Path] = None,
        workdir: Optional[Path] = None,
    ):
        """Create the workflow API.

        Note that if provided, this also changes to the provided workdir.
        It will change back to the previous working directory when the workflow API object is deleted.

        Arguments
        ---------
        config_settings: ConfigSettings -- The config settings for the workflow.
        resource_settings: ResourceSettings -- The resource settings for the workflow.
        storage_settings: StorageSettings -- The storage settings for the workflow.
        snakefile: Optional[Path] -- The path to the snakefile. If not provided, default locations will be tried.
        workdir: Optional[Path] -- The path to the working directory. If not provided, the current working directory will be used.
        """

        if config_settings is None:
            config_settings = ConfigSettings()
        if storage_settings is None:
            storage_settings = StorageSettings()
        if workflow_settings is None:
            workflow_settings = WorkflowSettings()
        if deployment_settings is None:
            deployment_settings = DeploymentSettings()
        if storage_provider_settings is None:
            storage_provider_settings = dict()

        self._check_is_in_context()

        self.setup_logger(mode=workflow_settings.exec_mode)

        self._check_default_storage_provider(storage_settings=storage_settings)

        snakefile = resolve_snakefile(snakefile)

        self._workflow_api = WorkflowApi(
            snakemake_api=self,
            snakefile=snakefile,
            workdir=workdir,
            config_settings=config_settings,
            resource_settings=resource_settings,
            storage_settings=storage_settings,
            workflow_settings=workflow_settings,
            deployment_settings=deployment_settings,
            storage_provider_settings=storage_provider_settings,
        )
        return self._workflow_api

    def _cleanup(self):
        """Cleanup the workflow."""
        if not self.output_settings.keep_logger:
            logger.cleanup()
        if self._workflow_api is not None:
            self._workflow_api._workdir_handler.change_back()
            if self._workflow_api._workflow_store is not None:
                self._workflow_api._workflow_store.tear_down()

    def deploy_sources(
        self,
        query: str,
        checksum: str,
        storage_settings: StorageSettings,
        storage_provider_settings: Dict[str, TaggedSettings],
    ):
        if (
            storage_settings.default_storage_provider is None
            or storage_settings.default_storage_prefix is None
        ):
            raise ApiError(
                "A default storage provider and prefix has to be set for deployment of "
                "sources."
            )

        self._check_default_storage_provider(storage_settings=storage_settings)

        plugin = StoragePluginRegistry().get_plugin(
            storage_settings.default_storage_provider
        )
        if not plugin.is_read_write():
            raise ApiError(
                f"Default storage provider {storage_settings.default_storage_provider} "
                "is not a read-write storage provider."
            )

        plugin_settings = storage_provider_settings.get(
            storage_settings.default_storage_provider
        ).get_settings(None)

        plugin.validate_settings(plugin_settings)

        provider_instance = plugin.storage_provider(
            local_prefix=storage_settings.local_storage_prefix,
            settings=plugin_settings,
            is_default=True,
        )
        query_validity = provider_instance.is_valid_query(query)
        if not query_validity:
            raise ApiError(
                f"Error when applying default storage provider "
                f"{storage_settings.default_storage_provider} to upload workflow "
                f"sources. {query_validity}"
            )
        storage_object = provider_instance.object(query)
        async_run(storage_object.managed_retrieve())
        with open(storage_object.local_path(), "rb") as f:
            obtained_checksum = hashlib.file_digest(f, "sha256").hexdigest()
        if obtained_checksum != checksum:
            raise ApiError(
                f"Checksum of retrieved sources ({obtained_checksum}) does not match "
                f"expected checksum ({checksum})."
            )
        with tarfile.open(storage_object.local_path(), "r") as tar:
            tar.extractall()

    def print_exception(self, ex: Exception):
        """Print an exception during workflow execution in a human readable way
        (with adjusted line numbers for exceptions raised in Snakefiles and stack
        traces that hide Snakemake internals for better readability).

        Arguments
        ---------
        ex: Exception -- The exception to print.
        """
        linemaps = dict()
        if (
            self._workflow_api is not None
            and self._workflow_api._workflow_store is not None
        ):
            linemaps = self._workflow_api._workflow_store.linemaps
        print_exception(ex, linemaps)

    def setup_logger(
        self,
        stdout: bool = False,
        mode: ExecMode = ExecMode.DEFAULT,
        dryrun: bool = False,
    ):
        if not self.output_settings.keep_logger:
            setup_logger(
                handler=self.output_settings.log_handlers,
                quiet=self.output_settings.quiet,
                nocolor=self.output_settings.nocolor,
                debug=self.output_settings.verbose,
                printshellcmds=self.output_settings.printshellcmds,
                debug_dag=self.output_settings.debug_dag,
                stdout=stdout or self.output_settings.stdout,
                mode=mode,
                show_failed_logs=self.output_settings.show_failed_logs,
                dryrun=dryrun,
            )

    def _check_is_in_context(self):
        if not self._is_in_context:
            raise ApiError(
                "This method can only be called when SnakemakeApi is used within a with "
                "statement."
            )

    def _check_default_storage_provider(self, storage_settings: StorageSettings):
        if storage_settings.default_storage_provider is not None:
            plugin = StoragePluginRegistry().get_plugin(
                storage_settings.default_storage_provider
            )
            if not plugin.is_read_write():
                raise ApiError(
                    f"Default storage provider {storage_settings.default_storage_provider} "
                    "is not a read-write storage provider."
                )

    def __enter__(self):
        self._is_in_context = True
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._is_in_context = False
        self._cleanup()


def _no_dag(method):
    def _handle_no_dag(self: "WorkflowApi", *args, **kwargs):
        self.resource_settings.cores = 1
        return method(self, *args, **kwargs)

    return _handle_no_dag


@dataclass
class WorkflowApi(ApiBase):
    """The workflow API.

    Arguments
    ---------
    snakemake_api: SnakemakeApi -- The Snakemake API.
    snakefile: Path -- The path to the snakefile.
    config_settings: ConfigSettings -- The config settings for the workflow.
    resource_settings: ResourceSettings -- The resource settings for the workflow.
    """

    snakemake_api: SnakemakeApi
    snakefile: Path
    workdir: Optional[Path]
    config_settings: ConfigSettings
    resource_settings: ResourceSettings
    storage_settings: StorageSettings
    workflow_settings: WorkflowSettings
    deployment_settings: DeploymentSettings
    storage_provider_settings: Mapping[str, TaggedSettings]

    _workflow_store: Optional[Workflow] = field(init=False, default=None)
    _workdir_handler: Optional[WorkdirHandler] = field(init=False)

    def dag(
        self,
        dag_settings: Optional[DAGSettings] = None,
    ):
        """Create a DAG API.

        Arguments
        ---------
        dag_settings: DAGSettings -- The DAG settings for the DAG API.
        """
        if dag_settings is None:
            dag_settings = DAGSettings()

        return DAGApi(
            self.snakemake_api,
            self,
            dag_settings=dag_settings,
        )

    @_no_dag
    def lint(self, json: bool = False):
        """Lint the workflow.

        Arguments
        ---------
        json: bool -- Whether to print the linting results as JSON.

        Returns
        -------
        True if any lints were printed
        """
        workflow = self._get_workflow(check_envvars=False)
        self._workflow_store = workflow
        workflow.include(
            self.snakefile, overwrite_default_target=True, print_compilation=False
        )
        workflow.check()
        return workflow.lint(json=json)

    @_no_dag
    def list_rules(self, only_targets: bool = False):
        """List the rules of the workflow.

        Arguments
        ---------
        only_targets: bool -- Whether to only list target rules.
        """
        self._workflow.list_rules(only_targets=only_targets)

    @_no_dag
    def list_resources(self):
        """List the resources of the workflow."""
        self._workflow.list_resources()

    @_no_dag
    def print_compilation(self):
        """Print the pure python compilation of the workflow."""
        workflow = self._get_workflow()
        workflow.include(self.snakefile, print_compilation=True)

    @property
    def _workflow(self):
        if self._workflow_store is None:
            workflow = self._get_workflow()
            self._workflow_store = workflow
            workflow.include(
                self.snakefile, overwrite_default_target=True, print_compilation=False
            )
            workflow.check()
        return self._workflow_store

    def _get_workflow(self, **kwargs):
        from snakemake.workflow import Workflow

        if "group_settings" not in kwargs:
            # just init with defaults, can be overwritten later
            kwargs["group_settings"] = GroupSettings()

        return Workflow(
            config_settings=self.config_settings,
            resource_settings=self.resource_settings,
            workflow_settings=self.workflow_settings,
            deployment_settings=self.deployment_settings,
            storage_settings=self.storage_settings,
            output_settings=self.snakemake_api.output_settings,
            overwrite_workdir=self.workdir,
            storage_provider_settings=self.storage_provider_settings,
            **kwargs,
        )

    def __post_init__(self):
        self._workdir_handler = None
        super().__post_init__()
        self.snakefile = self.snakefile.absolute()
        self._workdir_handler = WorkdirHandler(self.workdir)
        self._workdir_handler.change_to()

    def _check(self):
        if not self.snakefile.exists():
            raise ApiError(f'Snakefile "{self.snakefile}" not found.')


@dataclass
class DAGApi(ApiBase):
    """The DAG API.

    Arguments
    ---------
    snakemake_api: SnakemakeApi -- The Snakemake API.
    workflow_api: WorkflowApi -- The workflow API.
    dag_settings: DAGSettings -- The DAG settings for the DAG API.
    """

    snakemake_api: SnakemakeApi
    workflow_api: WorkflowApi
    dag_settings: DAGSettings

    def __post_init__(self):
        self.workflow_api._workflow.dag_settings = self.dag_settings

    def execute_workflow(
        self,
        executor: str = "local",
        execution_settings: Optional[ExecutionSettings] = None,
        remote_execution_settings: Optional[RemoteExecutionSettings] = None,
        scheduling_settings: Optional[SchedulingSettings] = None,
        group_settings: Optional[GroupSettings] = None,
        executor_settings: Optional[ExecutorSettingsBase] = None,
        updated_files: Optional[List[str]] = None,
    ):
        """Execute the workflow.

        Arguments
        ---------
        executor: str -- The executor to use.
        execution_settings: ExecutionSettings -- The execution settings for the workflow.
        resource_settings: ResourceSettings -- The resource settings for the workflow.
        remote_execution_settings: RemoteExecutionSettings -- The remote execution settings for the workflow.
        executor_settings: Optional[ExecutorSettingsBase] -- The executor settings for the workflow.
        updated_files: Optional[List[str]] -- An optional list where Snakemake will put all updated files.
        """

        if execution_settings is None:
            execution_settings = ExecutionSettings()
        if remote_execution_settings is None:
            remote_execution_settings = RemoteExecutionSettings()
        if scheduling_settings is None:
            scheduling_settings = SchedulingSettings()
        if group_settings is None:
            group_settings = GroupSettings()

        if (
            remote_execution_settings.immediate_submit
            and not self.workflow_api.storage_settings.notemp
        ):
            raise ApiError(
                "immediate_submit has to be combined with notemp (it does not support temp file handling)"
            )

        executor_plugin_registry = ExecutorPluginRegistry()
        executor_plugin = executor_plugin_registry.get_plugin(executor)

        if executor_settings is not None:
            executor_plugin.validate_settings(executor_settings)

        if executor_plugin.common_settings.implies_no_shared_fs:
            # no shared FS at all
            self.workflow_api.storage_settings.shared_fs_usage = frozenset()

        if (
            executor_plugin.common_settings.local_exec
            and not executor_plugin.common_settings.dryrun_exec
            and self.workflow_api.workflow_settings.exec_mode == ExecMode.DEFAULT
        ):
            logger.info("Assuming unrestricted shared filesystem usage.")
            self.workflow_api.storage_settings.shared_fs_usage = SharedFSUsage.all()
        if executor_plugin.common_settings.job_deploy_sources:
            remote_execution_settings.job_deploy_sources = True

        if (
            self.workflow_api.workflow_settings.exec_mode == ExecMode.DEFAULT
            and SharedFSUsage.INPUT_OUTPUT
            not in self.workflow_api.storage_settings.shared_fs_usage
            and (
                not self.workflow_api.storage_settings.default_storage_provider
                or self.workflow_api.storage_settings.default_storage_prefix is None
            )
            and executor_plugin.common_settings.can_transfer_local_files is False
        ):
            raise ApiError(
                "If no shared filesystem is assumed for input and output files, a "
                "default storage provider (--default-storage-provider) and "
                "default storage prefix (--default-storage-prefix) has to be set. "
                "See https://snakemake.github.io/snakemake-plugin-catalog for possible "
                "storage provider plugins."
            )
        if (
            executor_plugin.common_settings.local_exec
            and not executor_plugin.common_settings.dryrun_exec
            and self.workflow_api.workflow_settings.exec_mode == ExecMode.DEFAULT
            and self.workflow_api.storage_settings.shared_fs_usage
            != SharedFSUsage.all()
        ):
            raise ApiError(
                "For local execution, --shared-fs-usage has to be unrestricted."
            )

        self.snakemake_api.setup_logger(
            stdout=executor_plugin.common_settings.dryrun_exec,
            mode=self.workflow_api.workflow_settings.exec_mode,
            dryrun=executor_plugin.common_settings.dryrun_exec,
        )

        if executor_plugin.common_settings.local_exec:
            if (
                not executor_plugin.common_settings.dryrun_exec
                and not executor_plugin.common_settings.touch_exec
            ):
                if self.workflow_api.resource_settings.cores is None:
                    raise ApiError(
                        "cores have to be specified for local execution "
                        "(use --cores N with N being a number >= 1 or 'all')"
                    )
                # clean up all previously recorded jobids.
                shell.cleanup()
            else:
                # set cores if that is not done yet
                if self.workflow_api.resource_settings.cores is None:
                    self.workflow_api.resource_settings.cores = 1
            if (
                execution_settings.debug
                and self.workflow_api.resource_settings.cores > 1
            ):
                raise ApiError(
                    "debug mode cannot be used with multi-core execution, "
                    "please enforce a single core by setting --cores 1"
                )
        else:
            if self.workflow_api.resource_settings.nodes is None:
                raise ApiError(
                    "maximum number of parallel jobs/used nodes has to be specified for remote execution "
                    "(use --jobs N with N being a number >= 1)"
                )
            # non local execution
            if self.workflow_api.resource_settings.default_resources is None:
                # use full default resources if in cluster or cloud mode
                self.workflow_api.resource_settings.default_resources = (
                    DefaultResources(mode="full")
                )
            if execution_settings.edit_notebook is not None:
                raise ApiError(
                    "notebook edit mode is only allowed with local execution."
                )
            if execution_settings.debug:
                raise ApiError("debug mode cannot be used with non-local execution")

        if executor_plugin.common_settings.touch_exec:
            # no actual execution happening, hence we can omit any deployment
            self.workflow_api.deployment_settings.deployment_method = frozenset()

        execution_settings.use_threads = (
            execution_settings.use_threads
            or (os.name not in ["posix"])
            or not executor_plugin.common_settings.local_exec
        )

        logger.setup_logfile()

        workflow = self.workflow_api._workflow
        workflow.execution_settings = execution_settings
        workflow.remote_execution_settings = remote_execution_settings
        workflow.scheduling_settings = scheduling_settings
        workflow.group_settings = group_settings

        workflow.execute(
            executor_plugin=executor_plugin,
            executor_settings=executor_settings,
            updated_files=updated_files,
        )

    def _no_exec(method):
        def _handle_no_exec(self, *args, **kwargs):
            self.workflow_api.resource_settings.cores = 1
            return method(self, *args, **kwargs)

        return _handle_no_exec

    @_no_exec
    def generate_unit_tests(self, path: Path):
        """Generate unit tests for the workflow.

        Arguments
        ---------
        path: Path -- The path to store the unit tests.
        """
        self.workflow_api._workflow.generate_unit_tests(path=path)

    @_no_exec
    def containerize(self):
        """Containerize the workflow."""
        self.workflow_api._workflow.containerize()

    @_no_exec
    def create_report(
        self,
        reporter: str = "html",
        report_settings: Optional[ReportSettingsBase] = None,
    ):
        """Create a report for the workflow.

        Arguments
        ---------
        report: Path -- The path to the report.
        report_stylesheet: Optional[Path] -- The path to the report stylesheet.
        reporter: str -- report plugin to use (default: html)
        """

        report_plugin_registry = ReportPluginRegistry()
        report_plugin = report_plugin_registry.get_plugin(reporter)

        if report_settings is not None:
            report_plugin.validate_settings(report_settings)

        self.workflow_api._workflow.create_report(
            report_plugin=report_plugin,
            report_settings=report_settings,
        )

    @_no_exec
    def printdag(self):
        """Print the DAG of the workflow."""
        self.workflow_api._workflow.printdag()

    @_no_exec
    def printrulegraph(self):
        """Print the rule graph of the workflow."""
        self.workflow_api._workflow.printrulegraph()

    @_no_exec
    def printfilegraph(self):
        """Print the file graph of the workflow."""
        self.workflow_api._workflow.printfilegraph()

    @_no_exec
    def printd3dag(self):
        """Print the DAG of the workflow in D3.js compatible JSON."""
        self.workflow_api._workflow.printd3dag()

    @_no_exec
    def unlock(self):
        """Unlock the workflow."""
        self.workflow_api._workflow.unlock()

    @_no_exec
    def cleanup_metadata(self, paths: List[Path]):
        """Cleanup the metadata of the workflow."""
        self.workflow_api._workflow.cleanup_metadata(paths)

    @_no_exec
    def conda_cleanup_envs(self):
        """Cleanup the conda environments of the workflow."""
        self.workflow_api.deployment_settings.imply_deployment_method(
            DeploymentMethod.CONDA
        )
        self.workflow_api._workflow.conda_cleanup_envs()

    @_no_exec
    def conda_create_envs(self):
        """Only create the conda environments of the workflow."""
        self.workflow_api.deployment_settings.imply_deployment_method(
            DeploymentMethod.CONDA
        )
        self.workflow_api._workflow.conda_create_envs()

    @_no_exec
    def conda_list_envs(self):
        """List the conda environments of the workflow."""
        self.workflow_api.deployment_settings.imply_deployment_method(
            DeploymentMethod.CONDA
        )
        self.workflow_api._workflow.conda_list_envs()

    @_no_exec
    def cleanup_shadow(self):
        """Cleanup the shadow directories of the workflow."""
        self.workflow_api._workflow.cleanup_shadow()

    @_no_exec
    def container_cleanup_images(self):
        """Cleanup the container images of the workflow."""
        self.workflow_api.deployment_settings.imply_deployment_method(
            DeploymentMethod.APPTAINER
        )
        self.workflow_api._workflow.container_cleanup_images()

    @_no_exec
    def list_changes(self, change_type: ChangeType):
        """List the changes of the workflow.

        Arguments
        ---------
        change_type: ChangeType -- The type of changes to list.
        """
        self.workflow_api._workflow.list_changes(change_type=change_type)

    @_no_exec
    def list_untracked(self):
        """List the untracked files of the workflow."""
        self.workflow_api._workflow.list_untracked()

    @_no_exec
    def summary(self, detailed: bool = False):
        """Summarize the workflow.

        Arguments
        ---------
        detailed: bool -- Whether to print a detailed summary.
        """
        self.workflow_api._workflow.summary(detailed=detailed)

    @_no_exec
    def archive(self, path: Path):
        """Archive the workflow.

        Arguments
        ---------
        path: Path -- The path to the archive.
        """
        self.workflow_api._workflow.archive(path=path)

    @_no_exec
    def delete_output(self, only_temp: bool = False, dryrun: bool = False):
        """Delete the output of the workflow.

        Arguments
        ---------
        only_temp: bool -- Whether to only delete temporary output.
        dryrun: bool -- Whether to only dry-run the deletion.
        """
        self.workflow_api._workflow.delete_output(only_temp=only_temp, dryrun=dryrun)

    def export_to_cwl(self, path: Path):
        """Export the workflow to CWL.

        Arguments
        ---------
        path: Path -- The path to the CWL file.
        """
        self.workflow_api._workflow.export_to_cwl(path=path)
