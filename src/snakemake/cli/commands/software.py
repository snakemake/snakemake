import click
from snakemake.cli.common import SnakemakeContext
from snakemake.cli.operations import (
    run_conda_cleanup_envs,
    run_conda_create_envs,
    run_conda_list_envs,
    run_container_cleanup_images,
)


@click.group()
def software():
    """Manage software deployment environments."""


@software.group()
def conda():
    """Manage conda environments."""


@conda.command(name="list")
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def conda_list(snakefile, directory):
    """List all conda environments and their locations."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_conda_list_envs(dag_api)


@conda.command(name="cleanup")
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def conda_cleanup(snakefile, directory):
    """Remove unused conda environments."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_conda_cleanup_envs(dag_api)


@conda.command(name="create")
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def conda_create(snakefile, directory):
    """Create conda environments without running the workflow."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_conda_create_envs(dag_api)


@software.group()
def container():
    """Manage container images."""


@container.command(name="cleanup")
@click.option("--snakefile", "-s", type=click.Path(), default=None)
@click.option("--directory", "-d", type=click.Path(), default=None)
def container_cleanup(snakefile, directory):
    """Remove unused container images."""
    ctx = SnakemakeContext(snakefile=snakefile, workdir=directory)
    with ctx.dag() as dag_api:
        run_container_cleanup_images(dag_api)
