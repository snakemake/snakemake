"""
Tests for Snakemake’s API
"""

import time
from snakemake.api import snakemake
import asyncio
import tempfile
import os.path
from textwrap import dedent

from snakemake.scheduler import JobRateLimiter
from snakemake.settings.types import MaxJobsPerTimespan


def test_keep_logger():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "Snakefile")
        with open(path, "w") as f:
            print("rule:\n  output: 'result.txt'\n  shell: 'touch {output}'", file=f)
        snakemake(path, workdir=tmpdir, keep_logger=True)


def test_workflow_calling():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "Snakefile")
        with open(path, "w") as f:
            print(
                dedent(
                    """
                rule:
                    output: 'result.txt'
                    run:
                        with open(output[0], 'w') as f:
                            print("hello", file=f)
                """
                ),
                file=f,
            )
        workflow = Workflow(snakefile=snakefile, overwrite_workdir=tmpdir)


def test_run_script_directive():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "Snakefile")
        with open(path, "w") as f:
            print(
                dedent(
                    """
                rule:
                    output: 'result.txt'
                    run:
                        with open(output[0], 'w') as f:
                            print("hello", file=f)
                """
                ),
                file=f,
            )
        snakemake(path, workdir=tmpdir)


def test_run_script_directive_async():
    """Tests :func`snakemake.common.async_run`. The test ensures the ability to
    execute Snakemake API even if an asyncio event loop is already running.

    """
    import tracemalloc
    from snakemake.common import async_run

    tracemalloc.start()

    async def dummy_task():
        await asyncio.sleep(0.00001)

    async def main():
        async_run(dummy_task())
        test_run_script_directive()

    asyncio.run(main())


def test_dicts_in_config():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "Snakefile")
        with open(path, "w") as f:
            print(
                dedent(
                    """
                rule:
                    output: 'result.txt'
                    run:
                        with open(output[0], 'w') as f:
                            print("hello, this option " + config["this_option"] + "; this test dictionary " + config["test"]["this_dict"], file=f)
                """
                ),
                file=f,
            )
        snakemake(
            path,
            workdir=tmpdir,
            config={
                "this_option": "does_not_break",
                "test": {"this_dict": "shoult_not_either"},
            },
        )


def test_lockexception():
    from snakemake.persistence import Persistence
    from snakemake.exceptions import LockException

    persistence = Persistence()
    persistence.all_inputfiles = lambda: ["A.txt"]
    persistence.all_outputfiles = lambda: ["B.txt"]
    with persistence.lock():
        try:
            persistence.lock()
        except LockException as e:
            return True
        assert False


def test_job_rate_limiter():
    rate_limiter = JobRateLimiter(MaxJobsPerTimespan.parse_choice("10/1h"))

    rate_limiter.register_jobs(5)
    assert rate_limiter.get_free_jobs() == 5

    rate_limiter.register_jobs(5)
    assert rate_limiter.get_free_jobs() == 0

    rate_limiter = JobRateLimiter(MaxJobsPerTimespan.parse_choice("9/1s"))

    rate_limiter.register_jobs(4)
    assert rate_limiter.get_free_jobs() == 5

    time.sleep(1)
    assert rate_limiter.get_free_jobs() == 9
