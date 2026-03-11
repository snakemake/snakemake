import click
from snakemake.cli.common import SnakemakeContext
from snakemake.cli.operations import run_lint


@click.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def lint(snakefile, directory):
    """Check the workflow for common issues."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.workflow() as workflow_api:
        if not run_lint(workflow_api):
            raise SystemExit(1)
