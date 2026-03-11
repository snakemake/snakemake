import click
from copier import run_copy
from plumbum.commands.processes import ProcessExecutionError

# TEMPLATE_DIR = "src/snakemake/templates/project"

DEFAULT_WORKFLOW_URL = (
    "https://github.com/snakemake-workflows/snakemake-workflow-template"
)


@click.command()
@click.argument(
    "name",
)
@click.option(
    "--from",
    "repo",
    help="Import a template from GitHub or Gitlab.",
)
def new(name, repo):

    if repo:
        try:
            run_copy(repo, name)
        except ProcessExecutionError:
            print(f"Repository {repo} could not be found.")
            return
    else:
        # Create a project from a the default workflow template repo
        run_copy(DEFAULT_WORKFLOW_URL, name)


if __name__ == "__main__":
    new()
