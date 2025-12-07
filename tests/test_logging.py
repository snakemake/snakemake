import os
import shutil
import sys
import subprocess as sp
import logging
from collections import Counter
from pathlib import Path
import json

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


def test_log_events_dryrun(caplog: pytest.LogCaptureFixture):
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


def test_log_events(
    caplog: pytest.LogCaptureFixture, capfd: pytest.CaptureFixture[str]
):
    """Test LogEvent counts of records captured during workflow run."""

    with caplog.at_level(logging.INFO):
        run(dpath("logging/test_logfile"), check_results=False)

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
        },
    )

    captured = capfd.readouterr()
    stderr_output = captured.err
    expected_in_stderr = [
        "Building DAG of jobs",
        "Job stats:",
        "Finished job",
        "localrule all:",
    ]

    for expected_msg in expected_in_stderr:
        assert (
            expected_msg in stderr_output
        ), f"Expected '{expected_msg}' not found in stderr output"


def test_rule_failure(caplog: pytest.LogCaptureFixture):
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
            LogEvent.ERROR: 3,
            LogEvent.JOB_ERROR: 6,
        },
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
            # These are filtered out with the default filter
            LogEvent.WORKFLOW_STARTED: 1 if has_filter else 0,
            LogEvent.DEBUG_DAG: None if has_filter else 0,
            # Only emitted if requested by plugin
            LogEvent.RULEGRAPH: 1 if needs_rulegraph else 0,
        },
    )
