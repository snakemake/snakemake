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
