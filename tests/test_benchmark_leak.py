"""Regression tests for benchmark monitor thread leaks.

ScheduledPeriodicTimer re-armed a fresh DaemonTimer on every fire without checking the
stop flag, so a cancel() landing between a fire and _action() was undone and the monitor
thread lived for the rest of the process. This fires on normal job completion, so at scale
the leaked threads accumulate until the scheduler starves.

See: https://github.com/snakemake/snakemake/issues/XXXX
"""

from unittest.mock import MagicMock

import pytest

from snakemake.benchmark import BenchmarkTimer, ScheduledPeriodicTimer, benchmarked


class _Counter(ScheduledPeriodicTimer):
    """A ScheduledPeriodicTimer whose work() just counts, no psutil needed."""

    def __init__(self):
        super().__init__(interval=1)
        self.fires = 0

    def work(self):
        self.fires += 1


def _make(stopped, timer):
    t = _Counter()
    t._stopped = stopped
    t._timer = timer
    t.fires = 0
    return t


def test_action_is_noop_once_stopped():
    """_action() must not re-arm once _stopped — closes the cancel/reschedule race."""
    t = _make(stopped=True, timer=None)

    t._action()

    assert t.fires == 0, "work() ran after stop"
    assert t._timer is None, "_action re-armed a cancelled timer (the leak)"


def test_cancel_is_none_safe_and_sets_flag_first():
    """cancel() before start() must not raise and must record the stop."""
    t = _make(stopped=False, timer=None)

    t.cancel()  # _timer is None here — must not AttributeError

    assert t._stopped is True


def test_benchmarked_cancels_when_body_raises():
    """A raising body must still cancel the monitor (running_time set in finally)."""
    with pytest.raises(ValueError), benchmarked() as record:
        raise ValueError("boom")

    # running_time defaults to None; the finally in benchmarked() sets it after cancel().
    assert record.running_time is not None, "monitor not stopped on exception path"


def test_work_stops_monitor_once_process_exited():
    """BenchmarkTimer.work() stops the monitor when the observed process is gone, so it
    doesn't keep polling a dead PID until the context closes."""
    main = MagicMock()
    main.is_running.return_value = False
    timer = BenchmarkTimer.__new__(BenchmarkTimer)
    timer.main = main
    timer._stopped = False
    timer._update_record = lambda: None  # isolate the stop check from actual sampling

    timer.work()

    assert timer._stopped is True


def test_action_does_not_rearm_when_work_self_terminates():
    """If work() stops the monitor (observed process exited), _action must not re-arm."""
    t = _make(stopped=False, timer=None)

    def stopping_work():
        t.fires += 1
        t._stopped = True  # as BenchmarkTimer.work does when the process is gone

    t.work = stopping_work
    t._action()

    assert t.fires == 1
    assert t._timer is None, "_action re-armed after work() self-terminated"
