import click
from snakemake.cli.common import SnakemakeContext
from snakemake.cli.operations import (
    run_printd3dag,
    run_printdag,
    run_printfilegraph,
    run_printrulegraph,
)


@click.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
@click.option(
    "--rulegraph",
    is_flag=True,
    default=False,
    help="Print rule dependency graph instead of job DAG.",
)
@click.option(
    "--filegraph", is_flag=True, default=False, help="Print file dependency graph."
)
@click.option(
    "--d3",
    "d3dag",
    is_flag=True,
    default=False,
    help="Print DAG in D3.js-compatible JSON.",
)
def dagviz(snakefile, directory, rulegraph, filegraph, d3dag):
    """Print the workflow DAG for visualization."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        if rulegraph:
            run_printrulegraph(dag_api)
        elif filegraph:
            run_printfilegraph(dag_api)
        elif d3dag:
            run_printd3dag(dag_api)
        else:
            run_printdag(dag_api)
