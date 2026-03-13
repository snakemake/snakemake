import click
from snakemake.cli.common import SnakemakeContext
from snakemake.cli.operations import (
    run_list_changes,
    run_list_rules,
    run_list_untracked,
    run_print_compilation,
    run_summary,
)
from snakemake.settings.enums import ChangeType


@click.group()
def info():
    """Inspect workflow rules, changes, and output files."""


@info.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def rules(snakefile, directory):
    """List available rules."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.workflow() as workflow_api:
        run_list_rules(workflow_api)


@info.command(name="target-rules")
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def target_rules(snakefile, directory):
    """List available target rules."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.workflow() as workflow_api:
        run_list_rules(workflow_api, only_targets=True)


@info.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
@click.option("--detailed", is_flag=True, default=False, help="Include input files and shell commands.")
def summary(snakefile, directory, detailed):
    """Print a summary of all output files."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_summary(dag_api, detailed=detailed)


@info.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
@click.argument("change_type", type=click.Choice(["code", "input", "params"]))
def changes(snakefile, directory, change_type):
    """List output files with changed code, input, or params."""
    type_map = {
        "code": ChangeType.CODE,
        "input": ChangeType.INPUT,
        "params": ChangeType.PARAMS,
    }
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_list_changes(dag_api, type_map[change_type])


@info.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def untracked(snakefile, directory):
    """List files in the workdir not used by the workflow."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_list_untracked(dag_api)


@info.command()
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def compilation(snakefile, directory):
    """Print the compiled Python representation of the workflow."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.workflow() as workflow_api:
        run_print_compilation(workflow_api)