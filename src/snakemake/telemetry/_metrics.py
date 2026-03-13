"""Per-job system state metrics for telemetry.

Captures system "distress" metrics like load average, memory availability,
swap usage, and I/O wait percentage at job execution time.
"""

from __future__ import annotations

import os
import subprocess
import sys


def read_loadavg() -> float:
    """Return 1-minute load average."""
    try:
        return os.getloadavg()[0]
    except OSError:
        return 0.0


def read_mem_swap() -> tuple[int, int]:
    """Return (mem_avail_mb, swap_used_mb)."""
    mem_avail = 0
    swap_used = 0
    try:
        if sys.platform == "linux":
            swap_total = swap_free = 0
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        mem_avail = int(line.split()[1]) // 1024  # kB -> MB
                    elif line.startswith("SwapTotal:"):
                        swap_total = int(line.split()[1]) // 1024
                    elif line.startswith("SwapFree:"):
                        swap_free = int(line.split()[1]) // 1024
            swap_used = swap_total - swap_free
        elif sys.platform == "darwin":
            # Memory available: use vm_stat (pages × page_size)
            try:
                out = subprocess.check_output(["vm_stat"], text=True, timeout=2)
                page_size = 4096  # default macOS page size
                free_pages = inactive_pages = 0
                for line in out.splitlines():
                    if "page size of" in line:
                        page_size = int(line.split()[-2])
                    elif "Pages free:" in line:
                        free_pages = int(line.split()[-1].rstrip("."))
                    elif "Pages inactive:" in line:
                        inactive_pages = int(line.split()[-1].rstrip("."))
                mem_avail = (free_pages + inactive_pages) * page_size // (1024 * 1024)
            except (subprocess.SubprocessError, FileNotFoundError):
                pass
            # Swap used: sysctl vm.swapusage
            try:
                out = subprocess.check_output(
                    ["sysctl", "-n", "vm.swapusage"], text=True, timeout=2
                )
                # Format: "total = 2048.00M  used = 123.45M  free = 1924.55M"
                parts = out.split()
                for i, p in enumerate(parts):
                    if p == "used":
                        val = parts[i + 2].rstrip("M")
                        swap_used = int(float(val))
                        break
            except (
                subprocess.SubprocessError,
                ValueError,
                FileNotFoundError,
                IndexError,
            ):
                pass
    except Exception:
        pass
    return mem_avail, swap_used


def read_iowait_ticks() -> tuple[int, int]:
    """Return (iowait_ticks, total_ticks) from /proc/stat. Linux only."""
    if sys.platform != "linux":
        return 0, 0
    try:
        with open("/proc/stat") as f:
            for line in f:
                if line.startswith("cpu "):
                    fields = line.split()
                    # user, nice, system, idle, iowait, irq, softirq, steal
                    vals = [int(x) for x in fields[1:9]]
                    iowait = vals[4] if len(vals) > 4 else 0
                    total = sum(vals)
                    return iowait, total
    except (OSError, ValueError):
        pass
    return 0, 0


def compute_iowait_pct(start: tuple[int, int], end: tuple[int, int]) -> float:
    """Compute iowait percentage from before/after tick snapshots."""
    d_iowait = end[0] - start[0]
    d_total = end[1] - start[1]
    if d_total <= 0:
        return 0.0
    return round(100.0 * d_iowait / d_total, 2)


def capture_system_snapshot() -> dict:
    """Capture system state metrics at a point in time.

    Returns a dict with load_avg, mem_avail_mb, swap_used_mb, and
    iowait_start (raw ticks for later delta computation).
    """
    mem_avail, swap_used = read_mem_swap()
    return {
        "load_avg": round(read_loadavg(), 2),
        "mem_avail_mb": mem_avail,
        "swap_used_mb": swap_used,
        "iowait_start": read_iowait_ticks(),
    }
