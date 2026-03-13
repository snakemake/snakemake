from pathlib import Path

import click
from snakemake.cli.common import SnakemakeContext
from snakemake.cli.operations import (
    run_archive,
    run_containerize,
    run_generate_unit_tests,
)


@click.group()
def utils():
    """Miscellaneous workflow utilities."""


@utils.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def containerize(snakefile, directory):
    """Print a Dockerfile for the workflow."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_containerize(dag_api)


@utils.command(name="generate-unit-tests")
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
@click.argument("path", type=click.Path(), default=".tests/unit")
def generate_unit_tests(snakefile, directory, path):
    """Generate unit tests for each workflow rule."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_generate_unit_tests(dag_api, Path(path))


@utils.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
@click.argument("file", type=click.Path(writable=True))
def archive(snakefile, directory, file):
    """Archive the workflow into a tar file."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_archive(dag_api, Path(file))