import os
from pathlib import Path
import subprocess
import sys
import tempfile

from .common import dpath, get_expected_files, md5sum, serve_directory


def test_remote_snakefile_multiple_includes():
    source_dir = dpath("test_multiple_includes")
    expected_results = source_dir / "expected-results"
    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [str(repo_root), env.get("PYTHONPATH")])
    )

    with tempfile.TemporaryDirectory(prefix="snakemake-remote-snakefile-") as workdir:
        workdir = Path(workdir)

        with serve_directory(source_dir) as server_url:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "snakemake",
                    "--snakefile",
                    f"{server_url}/Snakefile",
                    "--cores",
                    "1",
                ],
                cwd=workdir,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )

        assert result.returncode == 0, result.stdout

        for relpath in get_expected_files(expected_results):
            output = workdir / relpath
            expected = expected_results / relpath
            assert output.exists(), f"Missing output {relpath}\n{result.stdout}"
            assert md5sum(output) == md5sum(expected)
