import click
from snakemake.cli.common import SnakemakeContext
from snakemake.cli.operations import run_unlock


@click.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def unlock(snakefile, directory):
    """Remove the workflow lock."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_unlock(dag_api)
