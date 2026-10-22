import json
import os
import socket
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import time
from pathlib import Path
from typing import Any, Callable

import psutil
import pytest

from snakemake import api
from snakemake.api import ApiError
from snakemake.settings import types as settings

sys.path.insert(0, os.path.dirname(__file__))

from .common import dpath, prepare_tmpdir, run

TIMEOUT = 60


class DapClient:
    """Minimal Debug Adapter Protocol client, standing in for an IDE."""

    def __init__(self, sock: socket.socket) -> None:
        self.file = sock.makefile("rwb")
        self.seq = 0

    def send(self, command: str, **arguments: Any) -> None:
        self.seq += 1
        message = {
            "seq": self.seq,
            "type": "request",
            "command": command,
            "arguments": arguments,
        }
        body = json.dumps(message).encode()
        self.file.write(f"Content-Length: {len(body)}\r\n\r\n".encode() + body)
        self.file.flush()

    def wait(self, matches: Callable[[dict], bool]) -> dict:
        while True:
            headers = dict(
                line.decode().split(": ", 1)
                for line in iter(self.file.readline, b"\r\n")
            )
            message = json.loads(self.file.read(int(headers["Content-Length"])))
            if matches(message):
                return message

    def request(self, command: str, **arguments: Any) -> dict:
        self.send(command, **arguments)
        return self.wait(lambda m: m.get("request_seq") == self.seq)


def connect(port: int) -> socket.socket:
    deadline = time.monotonic() + TIMEOUT
    while True:
        try:
            return socket.create_connection(("127.0.0.1", port), timeout=TIMEOUT)
        except ConnectionRefusedError:
            if time.monotonic() > deadline:
                raise
            time.sleep(0.05)


def stop_at_breakpoint(port: int, source: Path, line: int) -> dict:
    """Attach like an IDE, break at source:line, continue and return the top frame."""
    with connect(port) as sock:
        client = DapClient(sock)
        client.request("initialize", adapterID="debugpy", pathFormat="path")
        client.send("attach")
        client.wait(lambda m: m.get("event") == "initialized")
        client.request(
            "setBreakpoints", source={"path": str(source)}, breakpoints=[{"line": line}]
        )
        client.send("configurationDone")
        stopped = client.wait(lambda m: m.get("event") == "stopped")
        thread_id = stopped["body"]["threadId"]
        trace = client.request("stackTrace", threadId=thread_id, levels=1)
        client.request("continue", threadId=thread_id)
        return trace["body"]["stackFrames"][0]


def free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def wait_until(condition: Callable[[], bool]) -> None:
    deadline = time.monotonic() + TIMEOUT
    while not condition():
        assert time.monotonic() < deadline
        time.sleep(0.05)


def listens_on(process: psutil.Process, port: int) -> bool:
    return any(
        conn.status == psutil.CONN_LISTEN and conn.laddr.port == port
        for child in process.children(recursive=True)
        for conn in child.net_connections()
    )


@pytest.mark.parametrize("host, warned", [("127.0.0.1", False), ("0.0.0.0", True)])
def test_debug_scripts_stops_at_breakpoint(
    host: str,
    warned: bool,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """An attached debugger stops at a breakpoint set in the real script file."""
    port = free_port()
    monkeypatch.setenv("SNAKEMAKE_DEBUGPY_HOST", host)
    monkeypatch.setenv("SNAKEMAKE_DEBUGPY_PORT", str(port))
    path = dpath("test_debug_scripts")
    with prepare_tmpdir(path) as tmpdir:
        source = (Path(tmpdir) / "scripts" / "hello.py").resolve()
        with ThreadPoolExecutor(1) as executor:
            client = executor.submit(stop_at_breakpoint, port, source, 2)
            run(
                path,
                tmpdir=tmpdir,
                cleanup=False,
                shellcmd="snakemake --cores 1 --debug-scripts",
                timeout=TIMEOUT,
            )
            frame = client.result(TIMEOUT)

    assert (Path(frame["source"]["path"]), frame["line"]) == (source, 2)
    assert ("reachable from other machines" in capsys.readouterr().out) == warned


def test_debug_scripts_job_exits_when_snakemake_is_killed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A script waiting for a debugger exits once the Snakemake process is gone."""
    port = free_port()
    monkeypatch.setenv("SNAKEMAKE_DEBUGPY_PORT", str(port))
    with prepare_tmpdir(dpath("test_debug_scripts")) as tmpdir:
        with subprocess.Popen(
            [sys.executable, "-m", "snakemake", "--cores", "1", "--debug-scripts"],
            cwd=tmpdir,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ) as process:
            snakemake = psutil.Process(process.pid)
            wait_until(lambda: listens_on(snakemake, port))
            jobs = snakemake.children(recursive=True)
            snakemake.kill()
            process.wait()
            _, alive = psutil.wait_procs(jobs, timeout=10)

    assert alive == []
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", port))


def test_debug_does_not_pause_scripts() -> None:
    run(
        dpath("test_debug_scripts"),
        shellcmd="snakemake --cores 1 --debug",
        timeout=TIMEOUT,
    )


@pytest.mark.parametrize(
    "setting, option", [("debug", "--debug"), ("debug_scripts", "--debug-scripts")]
)
def test_debug_flags_reject_multicore_execution(setting: str, option: str) -> None:
    with api.SnakemakeApi(settings.OutputSettings()) as snakemake_api:
        dag_api = snakemake_api.workflow(
            resource_settings=settings.ResourceSettings(cores=2),
            snakefile=dpath("test_debug_scripts") / "Snakefile",
        ).dag()
        with pytest.raises(ApiError, match=f"^{option} cannot"):
            dag_api.execute_workflow(
                executor="local",
                execution_settings=settings.ExecutionSettings(**{setting: True}),
            )
