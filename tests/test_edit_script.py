import os
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path

import pytest
from jupyter_client.kernelspec import KernelSpecManager
from jupyter_client.manager import start_new_kernel

from snakemake import api
from snakemake.exceptions import WorkflowError
from snakemake.settings import types as settings

sys.path.insert(0, os.path.dirname(__file__))

from .common import dpath, prepare_tmpdir

TIMEOUT = 60


def wait_until(condition: Callable[[], bool]) -> None:
    deadline = time.monotonic() + TIMEOUT
    while not condition():
        assert time.monotonic() < deadline
        time.sleep(0.05)


def snakemake_kernels() -> list[str]:
    return [
        name
        for name in KernelSpecManager().find_kernel_specs()
        if name.startswith("snakemake-")
    ]


def test_edit_script_offers_kernel_then_runs_saved_script(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Editing a script offers a kernel with the job's snakemake object, then runs the saved script."""
    out = "out.txt"
    monkeypatch.setenv("JUPYTER_DATA_DIR", str(tmp_path))
    with prepare_tmpdir(dpath("test_edit_script")) as tmpdir:
        script = (Path(tmpdir) / "scripts" / "hello.py").resolve()
        process = subprocess.Popen(
            [sys.executable, "-m", "snakemake", "-c1", "--edit-script", out],
            cwd=tmpdir,
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        outputs = []
        try:
            wait_until(snakemake_kernels)
            [name] = snakemake_kernels()
            manager, client = start_new_kernel(kernel_name=name)
            client.execute_interactive(
                "import os; print(snakemake.output[0], os.getcwd(), __file__)",
                output_hook=lambda m: outputs.append(m["content"].get("text", "")),
                timeout=TIMEOUT,
            )
            manager.shutdown_kernel()
            script.write_text(script.read_text().replace('"hello"', '"edited"'))
            process.communicate("\n", timeout=TIMEOUT)
        finally:
            process.kill()
            process.wait()
        result = (Path(tmpdir) / out).read_text()

        assert "".join(outputs).split() == [out, os.path.realpath(tmpdir), str(script)]
    assert (process.returncode, result, snakemake_kernels()) == (0, "edited\n", [])


@pytest.mark.parametrize(
    ("target", "match"),
    [
        ("shell.txt", "requires a Python script"),
        ("remote.txt", "requires a local script"),
        ("container.txt", "does not support container jobs"),
    ],
)
def test_edit_script_rejects_rule_without_local_python_script(
    tmp_path: Path, target: str, match: str
) -> None:
    git = ["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@t"]
    (tmp_path / "hello.py").write_text("")
    subprocess.run([*git, "init"], check=True)
    subprocess.run([*git, "add", "hello.py"], check=True)
    subprocess.run([*git, "commit", "-m", "init"], check=True)
    with api.SnakemakeApi(settings.OutputSettings()) as snakemake_api:
        dag_api = snakemake_api.workflow(
            resource_settings=settings.ResourceSettings(cores=1),
            config_settings=settings.ConfigSettings(
                config={"remote_script": f"git+file:{tmp_path / 'hello.py'}@HEAD"}
            ),
            snakefile=dpath("test_edit_script") / "Snakefile",
        ).dag(settings.DAGSettings(targets={target}))
        with pytest.raises(WorkflowError, match=f"^--edit-script {match}"):
            dag_api.execute_workflow(
                executor="dryrun",
                execution_settings=settings.ExecutionSettings(edit_script=True),
            )
