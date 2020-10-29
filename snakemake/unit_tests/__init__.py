import textwrap
from itertools import groupby
from pathlib import Path
import shutil
import os

from snakemake.logging import logger
from snakemake import __version__


class RuleTest:
    def __init__(self, job, basedir):
        self.name = job.rule.name
        self.output = job.output
        self.path = basedir / self.name

    @property
    def target(self):
        return self.output or self.name

    @property
    def data_path(self):
        return self.path / "data"

    @property
    def expected_path(self):
        return self.path / "expected"


def generate(dag, path, deploy=["conda", "singularity"], configfiles=None):
    """Generate unit tests from given dag at given path."""
    logger.info("Generating unit tests for each rule...")

    try:
        from jinja2 import Template, Environment, PackageLoader
    except ImportError:
        raise WorkflowError(
            "Python package jinja2 must be installed to create reports."
        )

    env = Environment(
        loader=PackageLoader("snakemake", "unit_tests/templates"),
        trim_blocks=True,
        lstrip_blocks=True,
    )

    path = Path(path)
    os.makedirs(path, exist_ok=True)

    with open(path / "common.py", "w") as common:
        print(
            env.get_template("common.py.jinja2").render(version=__version__),
            file=common,
        )

    for rulename, jobs in groupby(dag.jobs, key=lambda job: job.rule.name):
        testpath = path / "test_{}.py".format(rulename)

        if testpath.exists():
            logger.info(
                "Skipping rule {} as a unit test already exists for it: {}.".format(
                    rulename, testpath
                )
            )
            continue

        written = False
        for job in jobs:
            if all(f.exists for f in job.input):
                logger.info(
                    "Generating unit test for rule {}: {}.".format(rulename, testpath)
                )
                os.makedirs(path / rulename, exist_ok=True)

                def copy_files(files, content_type):
                    for f in files:
                        f = Path(f)
                        target = path / rulename / content_type / f.parent
                        os.makedirs(target, exist_ok=True)
                        if f.is_dir():
                            shutil.copytree(f, target)
                        else:
                            shutil.copy(f, target)
                    if not files:
                        os.makedirs(path / rulename / content_type, exist_ok=True)
                        # touch gitempty file if there are no input files
                        open(path / rulename / content_type / ".gitempty", "w").close()

                copy_files(job.input, "data")
                copy_files(job.expanded_output, "expected")

                with open(testpath, "w") as test:
                    print(
                        env.get_template("ruletest.py.jinja2").render(
                            ruletest=RuleTest(job, path),
                            deploy=deploy,
                            configfiles=configfiles,
                        ),
                        file=test,
                    )

                written = True
                break

        if not written:
            logger.warning(
                "No job with all input files present for rule {}. "
                "Consider re-executing the workflow with --no-temp in order "
                "have all input files present before generating unit tests."
            )
