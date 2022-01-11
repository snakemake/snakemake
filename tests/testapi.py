"""
Tests for Snakemake’s API
"""
from snakemake import snakemake
import asyncio
import sys
import tempfile
import os.path
from textwrap import dedent


def test_keep_logger():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "Snakefile")
        with open(path, "w") as f:
            print("rule:\n  output: 'result.txt'\n  shell: 'touch {output}'", file=f)
        snakemake(path, workdir=tmpdir, keep_logger=True)


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

    if sys.version_info < (3, 7):
        async_run(main())
    else:
        asyncio.run(main())
