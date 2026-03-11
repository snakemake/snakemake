"""Workflow operations shared by all CLI frontends.

Each function performs a single workflow operation (lint, unlock, report, etc.)
given the appropriate API object. These are intentionally thin wrappers around
the Snakemake API — CLI-specific concerns (argument parsing, exit codes) belong
in the calling layer.
"""

from pathlib import Path
from typing import Optional

from snakemake.api import DAGApi, WorkflowApi
from snakemake.settings.enums import ChangeType
from snakemake.settings.types import GlobalReportSettings
from snakemake_interface_report_plugins.settings import ReportSettingsBase


def run_lint(workflow_api: WorkflowApi) -> bool:
    """Run linting. Returns True if the workflow is clean, False if issues found."""
    any_lint = workflow_api.lint()
    return not any_lint


def run_list_rules(workflow_api: WorkflowApi, only_targets: bool = False) -> None:
    """Print available rules (or only target rules) to stdout."""
    workflow_api.list_rules(only_targets=only_targets)


def run_print_compilation(workflow_api: WorkflowApi) -> None:
    """Print the compiled Python representation of the workflow."""
    workflow_api.print_compilation()


def run_unlock(dag_api: DAGApi) -> None:
    """Remove the workflow lock."""
    dag_api.unlock()


def run_printdag(dag_api: DAGApi) -> None:
    """Print the job DAG."""
    dag_api.printdag()


def run_printrulegraph(dag_api: DAGApi) -> None:
    """Print the rule dependency graph."""
    dag_api.printrulegraph()


def run_printfilegraph(dag_api: DAGApi) -> None:
    """Print the file dependency graph."""
    dag_api.printfilegraph()


def run_printd3dag(dag_api: DAGApi) -> None:
    """Print the DAG in D3.js-compatible JSON."""
    dag_api.printd3dag()


def run_summary(dag_api: DAGApi, detailed: bool = False) -> None:
    """Print a summary of all output files."""
    dag_api.summary(detailed=detailed)


def run_list_changes(dag_api: DAGApi, change_type: ChangeType) -> None:
    """List output files whose given property has changed."""
    dag_api.list_changes(change_type)


def run_list_untracked(dag_api: DAGApi) -> None:
    """List files in the workdir not tracked by the workflow."""
    dag_api.list_untracked()


def run_delete_output(
    dag_api: DAGApi, only_temp: bool = False, dryrun: bool = False
) -> None:
    """Delete workflow output files."""
    dag_api.delete_output(only_temp=only_temp, dryrun=dryrun)


def run_cleanup_metadata(dag_api: DAGApi, files: list[Path]) -> None:
    """Remove snakemake metadata for the given files."""
    dag_api.cleanup_metadata(files)


def run_cleanup_shadow(dag_api: DAGApi) -> None:
    """Remove stale shadow directories."""
    dag_api.cleanup_shadow()


def run_containerize(dag_api: DAGApi) -> None:
    """Print a Dockerfile for the workflow."""
    dag_api.containerize()


def run_generate_unit_tests(dag_api: DAGApi, path: Path) -> None:
    """Generate unit tests for each workflow rule."""
    dag_api.generate_unit_tests(path)


def run_archive(dag_api: DAGApi, path: Path) -> None:
    """Archive the workflow into a tar file."""
    dag_api.archive(path)


def run_report(
    dag_api: DAGApi,
    reporter: str,
    report_settings: ReportSettingsBase,
    metadata_template: Optional[Path] = None,
) -> None:
    """Generate a workflow report."""
    dag_api.create_report(
        reporter=reporter,
        report_settings=report_settings,
        global_report_settings=GlobalReportSettings(
            metadata_template=metadata_template,
        ),
    )


def run_conda_list_envs(dag_api: DAGApi) -> None:
    """List all conda environments and their locations."""
    dag_api.conda_list_envs()


def run_conda_cleanup_envs(dag_api: DAGApi) -> None:
    """Remove unused conda environments."""
    dag_api.conda_cleanup_envs()


def run_conda_create_envs(dag_api: DAGApi) -> None:
    """Create conda environments without running the workflow."""
    dag_api.conda_create_envs()


def run_container_cleanup_images(dag_api: DAGApi) -> None:
    """Remove unused container images."""
    dag_api.container_cleanup_images()
