from itertools import groupby
from pathlib import Path
import shutil
from snakemake.common import async_run

from snakemake.logging import logger
from snakemake import __version__
from snakemake.exceptions import WorkflowError


class RuleTest:
    def __init__(self, job, basedir):
        self.name = job.rule.name
        self.output = job.output
        self.path = basedir / self.name

    @property
    def target(self):
        return self.output or self.name

    @property
    def config_path(self):
        return self.path / "config"

    @property
    def data_path(self):
        return self.path / "data"

    @property
    def expected_path(self):
        return self.path / "expected"


def generate(
    dag, path: Path, deploy=None, snakefile=None, configfiles=None, rundir=None
):
    """Generate unit tests from given dag at a given path."""
    logger.info("Generating unit tests for each rule...")

    def copy_files(files, path):
        for f in set(files):
            f = Path(f)
            parent = f.parent
            if parent.is_absolute():
                root = str(f.parents[len(f.parents) - 1])
                parent = str(parent)[len(root) :]
            target = path / parent
            if f.is_dir():
                shutil.copytree(f, target / f.name)
            else:
                target.mkdir(parents=True, exist_ok=True)
                shutil.copy(f, target)
                (target / f.name).chmod(0o444)
        if not files:
            path.mkdir(parents=True, exist_ok=True)
            # touch gitempty file if there are no input files
            open(path / ".gitempty", "w").close()

    try:
        from jinja2 import Environment, PackageLoader
    except ImportError:
        raise WorkflowError(
            "Python package jinja2 must be installed to create reports."
        )

    env = Environment(
        loader=PackageLoader("snakemake", "unit_tests/templates"),
        trim_blocks=True,
        lstrip_blocks=True,
    )

    path.mkdir(parents=True, exist_ok=True)

    with open(path / "common.py", "w") as common:
        print(
            env.get_template("common.py.jinja2").render(version=__version__),
            file=common,
        )

    with open(path / "conftest.py", "w") as conftest:
        print(
            env.get_template("conftest.py.jinja2").render(version=__version__),
            file=conftest,
        )

    for rulename, jobs in groupby(dag.jobs, key=lambda job: job.rule.name):
        jobs = list(jobs)
        if jobs[0].rule.norun:
            logger.info(
                f"Skipping rule {rulename} because it does not execute anything."
            )
            continue

        testpath = path / f"test_{rulename}.py"

        if testpath.exists():
            logger.info(
                f"Skipping rule {rulename} as a unit test already exists for it: {testpath}."
            )
            continue

        written = False
        for job in jobs:
            if all(async_run(f.exists()) for f in job.input):
                logger.info(f"Generating unit test for rule {rulename}: {testpath}.")
                (path / rulename).mkdir(parents=True, exist_ok=True)

                copy_files(list(Path().glob("config*")), path / rulename / "config")
                copy_files(job.input, path / rulename / "data")
                copy_files(job.output, path / rulename / "expected")

                with open(testpath, "w") as test:
                    print(
                        env.get_template("ruletest.py.jinja2").render(
                            version=__version__,
                            ruletest=RuleTest(job, path.absolute().relative_to(rundir)),
                            deploy=deploy,
                            snakefile=Path(snakefile).absolute().relative_to(rundir),
                            configfiles=[
                                Path(config).absolute().relative_to(rundir)
                                for config in configfiles
                            ],
                        ),
                        file=test,
                    )

                written = True
                break

        if not written:
            logger.warning(
                "No job with all input files present for rule {}. "
                "Consider re-executing the workflow with --notemp in order "
                "have all input files present before generating unit tests."
            )
