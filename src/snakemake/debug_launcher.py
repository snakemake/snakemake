"""Debug launcher run in the job's Python; it must only import from the standard library."""

import os
import sys
import threading
import time
from functools import partial
from typing import Any, Callable, Dict


def _wait_posix(parent_pid: int) -> None:
    """Block until the parent process no longer exists (POSIX).

    Polls the parent with signal 0 every 0.5 seconds. A ``PermissionError``
    means the PID now belongs to another user's process, so the parent is gone.

    Args:
        parent_pid: PID of the Snakemake process that started the job.
    """
    while True:
        try:
            os.kill(parent_pid, 0)
        except (ProcessLookupError, PermissionError):
            return
        time.sleep(0.5)


def _windows_waiter(parent_pid: int) -> Callable[[], None]:
    """Create a function that blocks until the parent process exits (Windows).

    Opens the parent with ``SYNCHRONIZE`` access only. The open handle also
    prevents Windows from reusing the PID while the job runs.

    Args:
        parent_pid: PID of the Snakemake process that started the job.

    Returns:
        A function without arguments that waits on the parent process handle.

    Raises:
        OSError: If the parent process cannot be opened.
    """
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel32.WaitForSingleObject.restype = wintypes.DWORD
    kernel32.WaitForSingleObject.argtypes = (wintypes.HANDLE, wintypes.DWORD)
    handle = kernel32.OpenProcess(0x00100000, False, parent_pid)
    if not handle:
        raise ctypes.WinError(ctypes.get_last_error())
    return partial(kernel32.WaitForSingleObject, handle, 0xFFFFFFFF)


def _start_watchdog(parent_pid: int) -> None:
    """Start a daemon thread that ends the job when the parent process exits.

    This prevents a job that waits for a debugger from outliving a killed
    Snakemake process and keeping the debug port open.

    Args:
        parent_pid: PID of the Snakemake process that started the job.
    """
    wait = (
        _windows_waiter(parent_pid)
        if os.name == "nt"
        else partial(_wait_posix, parent_pid)
    )

    def watch() -> None:
        """Wait for the parent process to exit, then end this process."""
        wait()
        os._exit(1)

    threading.Thread(target=watch, daemon=True).start()


def launch(
    src: str,
    file_override: str,
    host: str,
    port: int,
    parent_pid: int,
    debugpy_dir: str,
    script_globals: Dict[str, Any],
) -> None:
    """Wait for an IDE debugger to attach, then run the job's script.

    The script is compiled with its real path as filename, so that breakpoints
    set in that file in the IDE are hit and line numbers match.

    Args:
        src: Path of the local script file to run, either the local script
            itself or its copy in the Snakemake source cache.
        file_override: Value of ``__file__`` seen by the script.
        host: Host on which debugpy listens.
        port: Port on which debugpy listens.
        parent_pid: PID of the Snakemake process that started the job.
        debugpy_dir: Directory containing the debugpy package of the Snakemake
            environment. It is appended to ``sys.path`` so that debugpy need not
            be installed in the job's environment.
        script_globals: Globals of the generated script, which contain the
            ``snakemake`` object defined by the preamble.
    """
    _start_watchdog(parent_pid)
    os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"
    sys.path.append(debugpy_dir)
    import debugpy

    debugpy.listen((host, port), in_process_debug_adapter=True)
    debugpy.wait_for_client()
    with open(src) as f:
        code = compile(f.read(), src, "exec")
    exec(code, {**script_globals, "__file__": file_override})
