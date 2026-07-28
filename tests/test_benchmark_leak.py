"""Regression tests for benchmark monitor thread leaks.

ScheduledPeriodicTimer re-armed a fresh DaemonTimer on every fire without checking the
stop flag, so a cancel() landing between a fire and _action() was undone and the monitor
thread lived for the rest of the process. This fires on normal job completion, so at scale
the leaked threads accumulate until the scheduler starves.
"""

from unittest.mock import MagicMock

import pytest

from snakemake import benchmark as bm
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


def test_action_does_not_rearm_when_cancelled_mid_action():
    """A cancel() that lands after _action's initial check (here: during work(), before
    the reschedule) must still prevent re-arming — the locked re-check catches it."""
    t = _make(stopped=False, timer=None)

    def work_then_cancel():
        t.fires += 1
        t.cancel()  # concurrent cancel racing the reschedule

    t.work = work_then_cancel
    t._action()

    assert t.fires == 1
    assert (
        t._timer is None
    ), "re-armed despite a cancel during the check→reschedule window"


def test_cancel_sets_flag_before_cancelling_timer():
    """cancel() must set _stopped before touching the timer, so a racing _action sees it."""
    t = _make(stopped=False, timer=None)
    observed = {}
    fake_timer = MagicMock()
    fake_timer.cancel.side_effect = lambda: observed.update(stopped=t._stopped)
    t._timer = fake_timer

    t.cancel()

    assert observed["stopped"] is True, "timer cancelled before the stop flag was set"
    fake_timer.cancel.assert_called_once()


def test_cancel_is_none_safe():
    """cancel() before start() (no timer yet) must not raise and must record the stop."""
    t = _make(stopped=False, timer=None)

    t.cancel()  # _timer is None here — must not AttributeError

    assert t._stopped is True


def test_benchmarked_cancels_monitor_even_when_body_raises(monkeypatch):
    """benchmarked() must cancel the monitor on the exception path (cancel in finally)."""
    fake = MagicMock()
    monkeypatch.setattr(bm, "BenchmarkTimer", lambda *a, **k: fake)

    with pytest.raises(ValueError), benchmarked():
        raise ValueError("boom")

    fake.cancel.assert_called_once()


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
