import os
import signal
import subprocess
import sys
import time
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest
from jupyter_client.kernelspec import KernelSpecManager
from jupyter_client.manager import start_new_kernel

from snakemake import api
from snakemake.exceptions import ApiError, WorkflowError
from snakemake.settings import types as settings

sys.path.insert(0, os.path.dirname(__file__))

from .common import dpath, prepare_tmpdir
from .conftest import skip_on_windows

TIMEOUT = 60
OUT = "out.txt"


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


@pytest.fixture
def paused_job(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> Iterator[tuple[subprocess.Popen, int, Path]]:
    import pty

    monkeypatch.setenv("JUPYTER_DATA_DIR", str(tmp_path))
    with prepare_tmpdir(dpath("test_edit_script")) as tmpdir:
        terminal, stdin = pty.openpty()
        process = subprocess.Popen(
            [sys.executable, "-m", "snakemake", "-c1", "--edit-script", OUT],
            cwd=tmpdir,
            stdin=stdin,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        os.close(stdin)
        try:
            wait_until(snakemake_kernels)
            yield process, terminal, Path(tmpdir)
        finally:
            process.kill()
            process.wait()
            os.close(terminal)


@skip_on_windows
def test_edit_script_offers_kernel_then_runs_saved_script(
    paused_job: tuple[subprocess.Popen, int, Path],
) -> None:
    """Editing a script offers a kernel with the job's snakemake object, then runs the saved script."""
    process, terminal, tmpdir = paused_job
    script = (tmpdir / "scripts" / "hello.py").resolve()
    [name] = snakemake_kernels()
    manager, client = start_new_kernel(kernel_name=name)
    outputs = []
    client.execute_interactive(
        "import os; print(snakemake.output[0], os.getcwd(), __file__)",
        output_hook=lambda m: outputs.append(m["content"].get("text", "")),
        timeout=TIMEOUT,
    )
    manager.shutdown_kernel()
    script.write_text(script.read_text().replace('"hello"', '"edited"'))
    os.write(terminal, b"\n")

    assert process.wait(timeout=TIMEOUT) == 0
    assert "".join(outputs).split() == [OUT, os.path.realpath(tmpdir), str(script)]
    assert ((tmpdir / OUT).read_text(), snakemake_kernels()) == ("edited\n", [])


@skip_on_windows
def test_edit_script_aborts_on_interrupt(
    paused_job: tuple[subprocess.Popen, int, Path],
) -> None:
    """Ctrl-C while paused fails the job without running the script and removes the kernel."""
    process, _, tmpdir = paused_job
    process.send_signal(signal.SIGINT)

    assert process.wait(timeout=TIMEOUT) != 0
    assert ((tmpdir / OUT).exists(), snakemake_kernels()) == (False, [])


@skip_on_windows
def test_unlock_removes_kernel_left_by_killed_job(
    paused_job: tuple[subprocess.Popen, int, Path],
) -> None:
    """A kernel left behind by a killed job is removed by --unlock."""
    process, _, tmpdir = paused_job
    process.kill()
    process.wait()
    left_behind = snakemake_kernels()
    subprocess.run(
        [sys.executable, "-m", "snakemake", "--unlock"], cwd=tmpdir, check=True
    )

    assert (len(left_behind), snakemake_kernels()) == (1, [])


@pytest.mark.parametrize(
    ("target", "executor", "cores", "error", "match"),
    [
        ("shell.txt", "dryrun", 1, WorkflowError, "requires a Python script"),
        ("remote.txt", "dryrun", 1, WorkflowError, "requires a local script"),
        ("container.txt", "dryrun", 1, WorkflowError, "does not support container"),
        (OUT, "dryrun", 2, ApiError, "cannot be used with multi-core execution"),
        (OUT, "local", 1, ApiError, "waits for Enter on stdin"),
    ],
)
def test_edit_script_rejects_invalid_setup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    target: str,
    executor: str,
    cores: int,
    error: type[Exception],
    match: str,
) -> None:
    monkeypatch.chdir(tmp_path)
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t"]
    (tmp_path / "hello.py").write_text("")
    subprocess.run([*git, "init"], check=True)
    subprocess.run([*git, "add", "hello.py"], check=True)
    subprocess.run([*git, "commit", "-m", "init"], check=True)
    with api.SnakemakeApi(settings.OutputSettings()) as snakemake_api:
        dag_api = snakemake_api.workflow(
            resource_settings=settings.ResourceSettings(cores=cores),
            config_settings=settings.ConfigSettings(
                config={"remote_script": f"git+file:{tmp_path / 'hello.py'}@HEAD"}
            ),
            snakefile=dpath("test_edit_script") / "Snakefile",
        ).dag(settings.DAGSettings(targets={target}))
        with pytest.raises(error, match=f"^--edit-script {match}"):
            dag_api.execute_workflow(
                executor=executor,
                execution_settings=settings.ExecutionSettings(edit_script=True),
            )
