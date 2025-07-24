import os
import shutil
import sys
import subprocess as sp
import pytest
import logging
from collections import Counter
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


def test_logging_config():
    """Relevant issue: https://github.com/snakemake/snakemake/issues/3044"""
    snakefile = os.path.join(dpath("logging/test_logging_config"), "Snakefile")
    p = sp.Popen(
        f"snakemake -s {snakefile}",
        shell=True,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
    )
    stdout, stderr = p.communicate()

    stdout = stdout.decode()
    assert p.returncode == 1
    assert "[TESTLOGGINGCONFIG]" in stdout


def test_logger_in_workflow():
    """relevant issue: https://github.com/snakemake/snakemake/issues/3558"""
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


def test_log_events_dryrun(caplog):
    with caplog.at_level(logging.INFO):
        run(dpath("logging/test_logfile"), executor="dryrun", check_results=False)

    events_found = []
    for record in caplog.records:
        if hasattr(record, "event") and record.event:
            events_found.append(record.event)

    event_counts = Counter(events_found)

    expected_event_counts = {
        LogEvent.WORKFLOW_STARTED: 1,
        LogEvent.RUN_INFO: 2,
        LogEvent.JOB_INFO: 6,
        LogEvent.SHELLCMD: 6,
    }

    for expected_event, expected_count in expected_event_counts.items():
        actual_count = event_counts.get(expected_event, 0)
        assert actual_count == expected_count, (
            f"Expected {expected_count} {expected_event} events, got {actual_count}. "
            f"All event counts: {event_counts}"
        )


def test_log_events(caplog, capfd):
    with caplog.at_level(logging.INFO):
        run(dpath("logging/test_logfile"), check_results=False)

    events_found = []
    for record in caplog.records:
        if hasattr(record, "event") and record.event:
            events_found.append(record.event)

    event_counts = Counter(events_found)

    expected_event_counts = {
        LogEvent.WORKFLOW_STARTED: 1,
        LogEvent.RUN_INFO: 1,
        LogEvent.JOB_INFO: 6,
        LogEvent.SHELLCMD: 6,
        LogEvent.RESOURCES_INFO: 2,
        LogEvent.PROGRESS: 6,
        LogEvent.JOB_STARTED: 4,
        LogEvent.JOB_FINISHED: 6,
    }

    # Check for unexpected events
    unexpected_events = set(event_counts.keys()) - set(expected_event_counts.keys())
    assert not unexpected_events, (
        f"Unexpected log events found: {unexpected_events}. "
        f"All event counts: {event_counts}"
    )

    # Check expected event counts
    for expected_event, expected_count in expected_event_counts.items():
        actual_count = event_counts.get(expected_event, 0)
        assert actual_count == expected_count, (
            f"Expected {expected_count} {expected_event} events, got {actual_count}. "
            f"All event counts: {event_counts}"
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


def test_rule_failure(caplog, capfd):
    run(dpath("logging/test_rule_failure"), check_results=False)
