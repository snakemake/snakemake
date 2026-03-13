from pathlib import Path

import click
from snakemake.cli.common import SnakemakeContext
from snakemake.cli.operations import (
    run_cleanup_metadata,
    run_cleanup_shadow,
    run_delete_output,
)


@click.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
@click.option("--all", "all_output", is_flag=True, default=False, help="Remove all output files.")
@click.option("--temp", "temp_output", is_flag=True, default=False, help="Remove only temporary output files.")
@click.option("--shadow", is_flag=True, default=False, help="Remove stale shadow directories.")
@click.option("--metadata", multiple=True, type=click.Path(), help="Remove metadata for the given files.")
@click.option("--dry-run", "-n", is_flag=True, default=False, help="List files without deleting.")
def clean(snakefile, directory, all_output, temp_output, shadow, metadata, dry_run):
    """Remove output files, shadow dirs, or metadata."""
    if not any([all_output, temp_output, shadow, metadata]):
        raise click.UsageError("Specify at least one of --all, --temp, --shadow, or --metadata.")

    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        if all_output:
            run_delete_output(dag_api, dryrun=dry_run)
        elif temp_output:
            run_delete_output(dag_api, only_temp=True, dryrun=dry_run)
        if shadow:
            run_cleanup_shadow(dag_api)
        if metadata:
            run_cleanup_metadata(dag_api, [Path(f) for f in metadata])