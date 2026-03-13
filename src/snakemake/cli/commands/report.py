from pathlib import Path

import click
from snakemake.cli.common import SnakemakeContext
from snakemake.cli.operations import run_report


@click.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
@click.option("--reporter", default="html", help="Report plugin to use.")
@click.option(
    "--metadata",
    type=click.Path(),
    default=None,
    help="Custom metadata file for the report.",
)
def report(snakefile, directory, reporter, metadata):
    """Generate a workflow report."""
    # TODO: once we figure out how to get plugin args
    pass
