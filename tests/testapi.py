"""
Tests for Snakemake’s API
"""
from snakemake import snakemake
import tempfile
import os.path
from textwrap import dedent


def test_keep_logger():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, 'Snakefile')
        with open(path, 'w') as f:
            print("rule:\n  output: 'result.txt'\n  shell: 'touch {output}'", file=f)
        snakemake(path, workdir=tmpdir, keep_logger=True)


def test_run_script_directive():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, 'Snakefile')
        with open(path, 'w') as f:
            print(dedent("""
                rule:
                    output: 'result.txt'
                    run:
                        with open(output[0], 'w') as f:
                            print("hello", file=f)
                """), file=f)
        snakemake(path, workdir=tmpdir)
