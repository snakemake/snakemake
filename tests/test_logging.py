"""Tests of logging system.


LogEvent types
--------------

Several tests check that specific LogEvent types have been generated, and verify that
all were formatted correctly with no errors. Together they should cover nearly the
entire set of events.

test_log_events (and others):

* WORKFLOW_STARTED
* RUN_INFO
* JOB_INFO
* SHELLCMD
* RESOURCES_INFO
* PROGRESS
* JOB_STARTED
* JOB_FINISHED
* DEBUG_DAG

test_rule_failure:

* JOB_ERROR

test_plugin (formatted output not checked):

* RULEGRAPH

Not currently in any tests:

* GROUP_INFO
* GROUP_ERROR
* ERROR
"""

from collections.abc import Iterable
import os
import shutil
import sys
import subprocess as sp
import logging
import logging.handlers
from collections import Counter
from pathlib import Path
from queue import Queue
import json
import re

import pytest

from snakemake_interface_logger_plugins.common import LogEvent

sys.path.insert(0, os.path.dirname(__file__))

from .common import run, dpath, apptainer, connected
from .conftest import (
    skip_on_windows,
    only_on_windows,
    ON_WINDOWS,
    needs_strace,
    ON_MACOS,
)

LOGGING_ERROR_SENTINEL = "--- Logging error ---"


def check_no_formatting_errors(stderr: str):
    """Assert that no log formatting errors occurred during the test.

    Python's logging.Handler.handleError() writes '--- Logging error ---'
    to stderr whenever a formatter raises an exception. This check ensures
    that all log events were formatted successfully.
    """
    assert (
        LOGGING_ERROR_SENTINEL not in stderr
    ), "Log formatting error detected in stderr output."


def count_events(caplog: pytest.LogCaptureFixture) -> dict[LogEvent, int]:
    """Count number of captured records for each event type."""

    counts = Counter()

    for record in caplog.records:
        event = getattr(record, "event", None)
        if event:
            counts[event] += 1

    return counts


def check_event_counts(
    observed: dict[LogEvent, int], expected: dict[LogEvent, int | None]
):
    """Check that captured event counts match expected values.

    Parameters
    ----------
    observed
        Observed event counts from ``count_events()``.
    expected
        Expected event counts. A value of None means any nonzero value.
    """

    unexpected = {
        event: count for event, count in observed.items() if event not in expected
    }

    try:
        assert not unexpected, f"Unexpected log events found: {unexpected}."

        for event, expected_count in expected.items():
            actual_count = observed.get(event, 0)
            if expected_count is None:
                assert actual_count > 0, f"Expected at least one {event} event."
            else:
                assert (
                    actual_count == expected_count
                ), f"Expected {expected_count} {event} events, got {actual_count}."

    except AssertionError:
        # Print all event counts to stderr for debugging if any checks fail.
        print("\nObserved event counts:", file=sys.stderr)
        for event in sorted(observed):
            print(f"  {event}: {observed[event]}", file=sys.stderr)

        raise


def check_log_contains(stderr: str, patterns: Iterable[str | re.Pattern]):
    """Check that the given strings or regex patterns are present in logging output."""

    for pattern in patterns:
        if isinstance(pattern, str):
            assert (
                pattern in stderr
            ), f"Expected substring {pattern!r} not found in stderr output"
        elif isinstance(pattern, re.Pattern):
            assert (
                pattern.search(stderr) is not None
            ), f"Expected pattern {pattern.pattern!r} not found in stderr output"
        else:
            raise TypeError(
                f"Elements of patterns must be str or re.Pattern, got {type(pattern)}"
            )


def test_logfile():
    import glob

    tmpdir = run(dpath("logging/test_logfile"), cleanup=False, check_results=False)
    finished_stmt = """
Finished jobid: 0 (Rule: all)
6 of 6 steps (100%) done"""

    log_dir = os.path.join(tmpdir, ".snakemake", "log")
    assert os.path.exists(log_dir), f"Log directory {log_dir} not found"

    log_files = glob.glob(os.path.join(log_dir, "*.snakemake.log"))
    assert log_files, "No log files found"

    log_files.sort(key=os.path.getmtime, reverse=True)
    latest_log = log_files[0]

    with open(latest_log, "r") as f:
        log_content = f.read()

    assert (
        finished_stmt.strip() in log_content.strip()
    ), f"Expected statement not found in log file. Log content: {log_content}"

    shutil.rmtree(tmpdir, ignore_errors=ON_WINDOWS)


def test_issue4063():
    # since it's difficult to replicate the exact conditions,
    # we'll do a unit test of the relevant logging code instead of a full workflow run.
    # the issue was that if the logger is configured with printshellcmds=True,
    # shell.py calls logger.info(None, extra={"event": LogEvent.SHELLCMD, "cmd": ...}),
    # which causes a crash because the default formatter expects the first argument to be a string message.
    from snakemake.logging import ColorizingTextHandler
    from snakemake.logging import DefaultFormatter
    from snakemake.logging import DefaultFilter

    handler = ColorizingTextHandler(stream=sys.stdout)
    formatter = DefaultFormatter(quiet=set(), show_failed_logs=False)
    handler.setFormatter(formatter)

    log_filter = DefaultFilter(
        quiet=set(), debug_dag=False, verbose=False, dryrun=False, printshellcmds=True
    )
    handler.addFilter(log_filter)

    # make the logging error a hard error
    def handle_logging_error(record):
        raise sys.exc_info()[1]

    handler.handleError = handle_logging_error

    test_logger = logging.getLogger("foo")
    test_logger.setLevel(logging.INFO)
    test_logger.addHandler(handler)

    test_logger.info(None, extra={"event": LogEvent.SHELLCMD, "cmd": "echo 'bar'"})
    test_logger.removeHandler(handler)


def test_formatting_error_produces_sentinel(capfd: pytest.CaptureFixture[str]):
    """Verify that a log formatting error produces the sentinel string in stderr.

    This validates that check_no_formatting_errors() would detect real formatting bugs.
    """
    from snakemake.logging import ColorizingTextHandler, DefaultFormatter

    handler = ColorizingTextHandler(stream=sys.stderr)
    formatter = DefaultFormatter()
    handler.setFormatter(formatter)

    test_logger = logging.getLogger("test_formatting_error_sentinel")
    test_logger.setLevel(logging.ERROR)
    test_logger.addHandler(handler)

    try:
        # Fire a JOB_ERROR event with a deliberately broken extra dict
        # (missing required keys like 'rule_name') to trigger a formatting error.
        test_logger.error(
            "bad record",
            extra=dict(event=LogEvent.JOB_ERROR),
        )
    finally:
        test_logger.removeHandler(handler)
        # test_logger remains globally registered in logging.Logger.manager.loggerDict, but Python
        # has no public API to remove it afterward. Shouldn't be a problem.

    stderr = capfd.readouterr().err
    assert LOGGING_ERROR_SENTINEL in stderr


def test_logging_config(tmp_path: Path):
    """Test configuring logging using ``logging.config.dictConfig()`` in the Snakefile.

    Relevant issue: https://github.com/snakemake/snakemake/issues/3044
    """
    snakefile = dpath("logging/test_logging_config") / "Snakefile"
    p = sp.Popen(
        f"snakemake -s {snakefile}",
        shell=True,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        cwd=tmp_path,
    )
    stdout, stderr = p.communicate()

    stdout = stdout.decode()
    assert p.returncode == 1
    assert "[TESTLOGGINGCONFIG]" in stdout


def test_logger_in_workflow():
    """Test adding a handler to the ``logger`` global in the Snakefile.

    relevant issue: https://github.com/snakemake/snakemake/issues/3558
    """
    import glob

    tmpdir = run(
        dpath("logging/test_workflow_logger"), cleanup=False, check_results=False
    )
    stmts = ["TESTINFO", "TESTWARN", "TESTERROR"]

    log_dir = os.path.join(tmpdir, ".snakemake", "log")
    assert os.path.exists(log_dir), f"Log directory {log_dir} not found"

    log_files = glob.glob(os.path.join(log_dir, "*.snakemake.log"))
    assert log_files, "No log files found"

    log_files.sort(key=os.path.getmtime, reverse=True)
    latest_log = log_files[0]

    with open(latest_log, "r") as f:
        log_content = f.read()

    custom_log = os.path.join(tmpdir, "mylog.txt")

    with open(custom_log, "r") as f:
        custom_log_content = f.read()

    for stmt in stmts:
        assert (
            stmt.strip() in log_content.strip()
        ), f"Expected statement {stmt} not found in log file. Log content: {log_content}"
        assert (
            stmt.strip() in custom_log_content.strip()
        ), f"Expected statement {stmt} not found in log file. Custom Log content: {custom_log_content}"

    shutil.rmtree(tmpdir, ignore_errors=ON_WINDOWS)


def test_log_events_dryrun(
    caplog: pytest.LogCaptureFixture,
    capfd: pytest.CaptureFixture[str],
):
    """Test LogEvent counts of records captured during dry run."""

    with caplog.at_level(logging.INFO):
        run(dpath("logging/test_logfile"), executor="dryrun", check_results=False)

    event_counts = count_events(caplog)

    check_event_counts(
        event_counts,
        {
            LogEvent.WORKFLOW_STARTED: 1,
            LogEvent.RUN_INFO: 2,
            LogEvent.JOB_INFO: 6,
            LogEvent.SHELLCMD: 6,
        },
    )

    # Check no errors occurred in formatting records
    stderr = capfd.readouterr().err
    check_no_formatting_errors(stderr)


def test_log_events(
    caplog: pytest.LogCaptureFixture,
    capfd: pytest.CaptureFixture[str],
):
    """Test LogEvent counts of records captured during workflow run."""

    with caplog.at_level(logging.DEBUG):
        # printshellcmds and debug_dag are enabled so that the SHELLCMD and DEBUG_DAG
        # events are also visible in stderr (otherwise the default filter would drop
        # them), letting check_log_contains() below cover every event type checked here.
        run(
            dpath("logging/test_logfile"),
            check_results=False,
            printshellcmds=True,
            debug_dag=True,
        )

    event_counts = count_events(caplog)

    check_event_counts(
        event_counts,
        {
            LogEvent.WORKFLOW_STARTED: 1,
            LogEvent.RUN_INFO: 1,
            LogEvent.JOB_INFO: 6,
            LogEvent.SHELLCMD: 6,
            LogEvent.RESOURCES_INFO: 2,
            LogEvent.PROGRESS: 6,
            LogEvent.JOB_STARTED: None,
            LogEvent.JOB_FINISHED: 6,
            LogEvent.DEBUG_DAG: None,
        },
    )

    stderr = capfd.readouterr().err

    # Check no errors occurred in formatting records
    check_no_formatting_errors(stderr)

    # Check expected lines in log output
    check_log_contains(
        stderr,
        [
            "Building DAG of jobs",
            "SNAKEMAKE",  # WORKFLOW_STARTED
            "Job stats:",  # RUN_INFO
            "localrule all:",  # JOB_INFO
            "Shell command:",  # SHELLCMD
            "Provided cores:",  # RESOURCES_INFO
            re.compile(r"\d+ of \d+ steps \(\d+%\) done"),  # PROGRESS
            re.compile(r"Execute \d+ jobs"),  # JOB_STARTED
            "Finished jobid:",  # JOB_FINISHED
            "candidate job",  # DEBUG_DAG
        ],
    )


def test_rule_failure(
    caplog: pytest.LogCaptureFixture,
    capfd: pytest.CaptureFixture[str],
):
    """Test LogEvent counts of records captured during workflow run with a failing rule."""

    with caplog.at_level(logging.INFO):
        run(dpath("logging/test_rule_failure"), check_results=False, shouldfail=True)

    event_counts = count_events(caplog)

    check_event_counts(
        event_counts,
        {
            LogEvent.WORKFLOW_STARTED: 1,
            LogEvent.RUN_INFO: 1,
            LogEvent.JOB_INFO: 3,
            LogEvent.SHELLCMD: 3,
            LogEvent.RESOURCES_INFO: 2,
            LogEvent.JOB_STARTED: None,
            LogEvent.JOB_ERROR: 6,
        },
    )

    # Check no errors occurred in formatting records
    stderr = capfd.readouterr().err
    check_no_formatting_errors(stderr)

    # Check expected lines in log output
    check_log_contains(
        stderr,
        [
            "Error in rule a:",  # JOB_ERROR
        ],
    )


@pytest.mark.parametrize(
    "stream,has_formatter,has_filter,needs_rulegraph",
    [
        (False, False, False, False),
        (True, False, False, False),
        (False, True, False, False),
        (False, False, True, False),
        (False, False, False, True),
    ],
)
def test_plugin(
    stream: bool,
    has_formatter: bool,
    has_filter: bool,
    needs_rulegraph: bool,
    tmp_path: Path,
):
    """Test using a logger plugin.

    Adds the logging/plugins/ directory to PYTHONPATH so the plugin registry can detect the
    "snakemake_logger_plugin_test" package. The plugin outputs each record in JSON format on a
    single line, including the "event" attribute so we can check event counts as in the other tests.
    The first line is a special record/event that reports information about how Snakemake has
    configured the handler (such as whether the default formatter or filter were attached).

    Parameters
    ----------
    stream
        If True output to stream, otherwise to file.
    has_formatter
        Value plugin handler should return for the "has_formatter" property.
    has_filter
        Value plugin handler should return for the "has_filter" property.
    needs_rulegraph
        Value plugin handler should return for the "needs_rulegraph" property.
    """

    plugin_dir = dpath("logging/plugins")
    test_dir = dpath("logging/test_logfile")
    outfile = tmp_path / "out.log"

    # Add plugin directory to PYTHONPATH so the registry can discover it
    env = dict(os.environ)
    current_path = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        str(plugin_dir)
        if current_path is None
        else str(plugin_dir) + os.pathsep + current_path
    )

    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "-s",
        str(test_dir / "Snakefile"),
        "-j1",
        "--verbose",
        "--printshellcmds",
        "--logger",
        "test",
    ]

    if not stream:
        cmd.extend(["--logger-test-outfile", outfile.name])
    if has_formatter:
        cmd.append("--logger-test-has-formatter")
    if has_filter:
        cmd.append("--logger-test-has-filter")
    if needs_rulegraph:
        cmd.append("--logger-test-needs-rulegraph")

    result = sp.run(
        cmd,
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=stream,
        text=True,
    )

    # Parse output file or stderr
    # If outputting to stderr it should replace Snakemake's default output.
    if stream:
        records = list(map(json.loads, result.stderr.splitlines()))
    else:
        with open(outfile) as fh:
            records = list(map(json.loads, fh))

    # Check logger info
    assert records[0]["event"] == "logger_info"
    assert records[0]["formatter_set"] == (not has_formatter)
    assert records[0]["filter_added"] == (not has_filter)

    # Check event counts
    event_counts = Counter(
        LogEvent[record["event"].upper()] for record in records[1:] if record["event"]
    )

    check_event_counts(
        event_counts,
        {
            LogEvent.RUN_INFO: 1,
            LogEvent.JOB_INFO: 6,
            LogEvent.SHELLCMD: 6,
            LogEvent.RESOURCES_INFO: 2,
            LogEvent.PROGRESS: 6,
            LogEvent.JOB_STARTED: None,
            LogEvent.JOB_FINISHED: 6,
            LogEvent.WORKFLOW_STARTED: 1,
            # Filtered out with the default filter
            LogEvent.DEBUG_DAG: None if has_filter else 0,
            # Only emitted if requested by plugin
            LogEvent.RULEGRAPH: 1 if needs_rulegraph else 0,
        },
    )


def test_stop_closes_handlers():
    """Test that LoggerManager.stop() closes and flushes all logging handlers.

    This includes plugin handlers managed by the QueueListener as well as those attached to the
    global logger instance.

    Regression test for https://github.com/snakemake/snakemake/issues/4136.
    Plugin handlers are attached to the QueueListener, not the logger itself, so they
    must be explicitly closed after the queue is drained.
    """
    from snakemake.logging import LoggerManager
    from snakemake.settings.types import OutputSettings

    class TrackingHandler(logging.Handler):
        """A handler that tracks whether emit(), flush(), and close() were called."""

        def __init__(self):
            super().__init__()
            self.records: list[logging.LogRecord] = []
            self.flush_called = False
            self.close_called = False

        def emit(self, record: logging.LogRecord) -> None:
            self.records.append(record)

        def flush(self) -> None:
            self.flush_called = True
            super().flush()

        def close(self) -> None:
            self.close_called = True
            super().close()

    test_logger = logging.getLogger("test_stop_closes_plugin_handlers")
    settings = OutputSettings(skip_plugin_handlers=True)
    manager = LoggerManager(test_logger, settings)

    # Simulate a plugin handler by manually wiring up a QueueListener.
    # This mirrors what _setup_plugins() does in production.
    plugin_handler = TrackingHandler()
    queue: Queue[logging.LogRecord] = Queue(-1)
    manager.queue_listener = logging.handlers.QueueListener(
        queue, plugin_handler, respect_handler_level=True
    )
    manager.queue_listener.start()
    manager.plugin_handlers = [plugin_handler]
    test_logger.addHandler(logging.handlers.QueueHandler(queue))

    # Also add a handler directly to the logger instance, similar to the default stream handler.
    global_handler = TrackingHandler()
    test_logger.addHandler(global_handler)

    try:
        test_logger.info("test message")
        manager.stop()
    finally:
        # Ensure the listener thread is stopped even if an assertion fails, to avoid
        # leaking the background thread into other tests.
        if manager.queue_listener._thread is not None:
            manager.queue_listener.stop()

    for handler, desc in [(plugin_handler, "plugin"), (global_handler, "global")]:
        assert (
            handler.flush_called
        ), f"{desc} handler was not flushed by LoggerManager.stop()"
        assert (
            handler.close_called
        ), f"{desc} handler was not closed by LoggerManager.stop()"
        assert (
            len(handler.records) == 1
        ), f"Expected 1 record delivered to {desc} handler, got {len(handler.records)}"
