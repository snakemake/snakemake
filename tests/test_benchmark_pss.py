"""Regression test for benchmark on Windows: missing psutil meminfo.pss.

On Windows, psutil.Process.memory_full_info() does not include the 'pss'
field (Proportional Set Size is Linux-only). BenchmarkTimer._update_record()
accessed meminfo.pss unconditionally, raising AttributeError which was
silently caught by work(), causing all benchmark metrics to be NA.

See: https://github.com/snakemake/snakemake/issues/4159
"""

from unittest.mock import MagicMock
import pytest
from snakemake.benchmark import BenchmarkTimer, BenchmarkRecord


def _make_mock_proc(has_pss):
    """Return a mock psutil.Process with or without pss in meminfo."""
    proc = MagicMock()
    proc.pid = 12345
    proc.name.return_value = "python.exe"

    if has_pss:
        meminfo = MagicMock(
            rss=100 * 1024 * 1024,
            vms=200 * 1024 * 1024,
            uss=80 * 1024 * 1024,
            pss=90 * 1024 * 1024,
        )
    else:
        meminfo = MagicMock(spec=["rss", "vms", "uss"])
        meminfo.rss = 100 * 1024 * 1024
        meminfo.vms = 200 * 1024 * 1024
        meminfo.uss = 80 * 1024 * 1024

    proc.memory_full_info.return_value = meminfo
    proc.cpu_percent.return_value = 50.0
    proc.cpu_times.return_value = MagicMock(user=1.0, system=0.5)
    proc.io_counters.return_value = MagicMock(
        read_bytes=1024, write_bytes=512
    )
    proc.children.return_value = []
    return proc


def test_update_record_without_pss():
    """Metrics must be collected even when meminfo lacks 'pss' (Windows)."""
    proc = _make_mock_proc(has_pss=False)

    record = BenchmarkRecord()
    timer = BenchmarkTimer.__new__(BenchmarkTimer)
    timer.pid = 12345
    timer.main = proc
    timer.bench_record = record
    timer.procs = {}

    timer._update_record()

    assert record.max_rss is not None, "max_rss is None — metrics not collected"
    assert record.max_vms is not None, "max_vms is None"
    assert record.max_uss is not None, "max_uss is None"
    assert record.max_rss == pytest.approx(100.0, abs=0.1)


def test_update_record_with_pss():
    """Sanity check: Linux-like meminfo with pss still works after fix."""
    proc = _make_mock_proc(has_pss=True)

    record = BenchmarkRecord()
    timer = BenchmarkTimer.__new__(BenchmarkTimer)
    timer.pid = 12345
    timer.main = proc
    timer.bench_record = record
    timer.procs = {}

    timer._update_record()

    assert record.max_rss == pytest.approx(100.0, abs=0.1)
    assert record.max_pss == pytest.approx(90.0, abs=0.1)