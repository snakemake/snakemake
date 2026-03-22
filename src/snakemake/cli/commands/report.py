"""Click subcommand for generating a Snakemake workflow report."""

from pathlib import Path
from typing import Optional

import click

from snakemake.cli.common import SnakemakeContext, workflow_options
from snakemake.cli.operations import run_report
from snakemake.cli.plugin_adapter import (
    add_plugin_options,
    plugin_settings_from_kwargs,
)
from snakemake_interface_report_plugins.registry import ReportPluginRegistry


@click.command()
@workflow_options
@click.option(
    "--reporter",
    type=click.Choice(list(ReportPluginRegistry().plugins.keys())),
    default="html",
    show_default=True,
    help="Report plugin to use.",
)
@click.option(
    "--report-metadata",
    type=click.Path(path_type=Path),
    default=None,
    help="Custom metadata file for the report landing page.",
)
@add_plugin_options("report")
def report(
    snakefile: Optional[Path],
    directory: Optional[Path],
    reporter: str,
    report_metadata: Optional[Path],
    **plugin_kwargs,
):
    """Generate a workflow report."""
    report_plugin = ReportPluginRegistry().get_plugin(reporter)
    report_settings = plugin_settings_from_kwargs(report_plugin, plugin_kwargs)

    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_report(
            dag_api,
            reporter=reporter,
            report_settings=report_settings,
            metadata_template=report_metadata,
        )