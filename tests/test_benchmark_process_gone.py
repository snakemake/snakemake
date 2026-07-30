"""A benchmarked process that exits before it can be sampled (a job faster than the
benchmark interval) is a benign, expected reason for missing stats — not a collection
failure. It is tracked separately from real errors and must not emit a warning.
"""

from unittest.mock import MagicMock

import psutil

from snakemake.benchmark import BenchmarkRecord, BenchmarkTimer


def test_update_record_marks_process_gone_on_nosuchprocess():
    """A vanished process is recorded as benign (process_gone), not appended to errors."""
    proc = MagicMock()
    proc.children.side_effect = psutil.NoSuchProcess(pid=999)

    record = BenchmarkRecord()
    timer = BenchmarkTimer.__new__(BenchmarkTimer)
    timer.main = proc
    timer.bench_record = record
    timer.procs = {}

    timer._update_record()

    assert record.process_gone is True
    assert record.errors == []


def _record_with(process_gone, errors):
    record = BenchmarkRecord()
    record.running_time = 0.01
    record.data_collected = False
    record.process_gone = process_gone
    record.errors = errors
    return record


def test_get_benchmarks_does_not_warn_when_only_reason_is_process_gone():
    record = _record_with(process_gone=True, errors=[])

    logger = MagicMock()
    import snakemake.benchmark as bm

    orig = bm.logger
    bm.logger = logger
    try:
        record.get_benchmarks()
    finally:
        bm.logger = orig

    logger.warning.assert_not_called()
    logger.debug.assert_called()  # downgraded to debug


def test_get_benchmarks_still_warns_on_real_collection_error():
    record = _record_with(process_gone=False, errors=["some real psutil/OS error"])

    logger = MagicMock()
    import snakemake.benchmark as bm

    orig = bm.logger
    bm.logger = logger
    try:
        record.get_benchmarks()
    finally:
        bm.logger = orig

    logger.warning.assert_called()  # genuine failures still surface


def _sampleable_proc(pid):
    """A mock process that yields one clean sample."""
    from types import SimpleNamespace

    proc = MagicMock()
    proc.pid = pid
    proc.memory_full_info.return_value = SimpleNamespace(
        rss=1024 * 1024, vms=1024 * 1024, uss=1024 * 1024
    )
    proc.io_counters.return_value = SimpleNamespace(read_bytes=0, write_bytes=0)
    proc.cpu_times.return_value = SimpleNamespace(user=0.1, system=0.1)
    proc.cpu_percent.return_value = 0.0
    proc.name.return_value = f"proc-{pid}"
    return proc


def _timer_for(main, children):
    main.children.return_value = children
    timer = BenchmarkTimer.__new__(BenchmarkTimer)
    timer.main = main
    timer.bench_record = BenchmarkRecord()
    timer.procs = {}
    return timer


def test_vanished_child_is_skipped_and_partial_sample_kept():
    """A child exiting mid-sample is normal churn: not process_gone, not an error,
    and the main process's sample still lands."""
    main = _sampleable_proc(999)
    child = _sampleable_proc(1000)
    child.oneshot.side_effect = psutil.NoSuchProcess(pid=1000)
    timer = _timer_for(main, [child])

    timer._update_record()

    assert timer.bench_record.process_gone is False
    assert timer.bench_record.errors == []
    assert timer.bench_record.data_collected is True
    assert timer.bench_record.max_rss == 1.0  # main's 1 MiB sampled despite the child


def test_main_vanishing_during_sampling_sets_process_gone():
    """The observed process disappearing mid-walk is the benign whole-job case."""
    main = _sampleable_proc(999)
    main.oneshot.side_effect = psutil.NoSuchProcess(pid=999)
    timer = _timer_for(main, [])

    timer._update_record()

    assert timer.bench_record.process_gone is True
    assert timer.bench_record.errors == []
    assert timer.bench_record.data_collected is False
