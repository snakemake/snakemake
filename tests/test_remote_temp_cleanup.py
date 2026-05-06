from collections import defaultdict
from unittest.mock import MagicMock

import pytest

# Import api first to resolve the dag <-> workflow circular import.
import snakemake.api  # noqa: F401
from snakemake.dag import DAG
from snakemake.io import flag


def _make_tempfile(name):
    return flag(name, "temp")


@pytest.fixture()
def mock_dag():
    """Create a minimal mock DAG with remote_exec=True to test is_needed_tempfile."""
    dag = object.__new__(DAG)

    dag.workflow = MagicMock()
    dag.workflow.subprocess_exec = False
    dag.workflow.remote_exec = True
    dag.workflow.storage_settings = MagicMock()
    dag.workflow.storage_settings.unneeded_temp_files = frozenset()

    dag._finished = set()
    dag._needrun = set()
    dag.depending = defaultdict(lambda: defaultdict(set))
    dag._derived_targetfiles = set()

    return dag


def test_temp_file_unneeded_when_no_downstream(mock_dag):
    producer = MagicMock()
    tempfile = _make_tempfile("output.tmp")

    mock_dag.depending[producer] = {}

    assert not mock_dag.is_needed_tempfile(producer, tempfile)


def test_temp_file_needed_when_downstream_unfinished(mock_dag):
    producer = MagicMock()
    consumer = MagicMock()
    tempfile = _make_tempfile("output.tmp")

    mock_dag.depending[producer] = {consumer: {tempfile}}
    mock_dag._needrun.add(consumer)

    assert mock_dag.is_needed_tempfile(producer, tempfile)


def test_temp_file_unneeded_when_downstream_finished(mock_dag):
    producer = MagicMock()
    consumer = MagicMock()
    tempfile = _make_tempfile("output.tmp")

    mock_dag.depending[producer] = {consumer: {tempfile}}
    mock_dag._finished.add(consumer)

    assert not mock_dag.is_needed_tempfile(producer, tempfile)


def test_temp_file_unneeded_remote_exec_empty_unneeded_set(mock_dag):
    """Regression: remote_exec with empty unneeded_temp_files must NOT prevent cleanup.

    Previously, is_needed_tempfile gated on unneeded_temp_files (always empty in the
    main process), causing temp files to never be cleaned from remote storage.
    """
    producer = MagicMock()
    tempfile = _make_tempfile("output.tmp")

    assert mock_dag.workflow.storage_settings.unneeded_temp_files == frozenset()
    mock_dag.depending[producer] = {}

    assert not mock_dag.is_needed_tempfile(producer, tempfile)


def test_subprocess_exec_always_returns_needed(mock_dag):
    mock_dag.workflow.subprocess_exec = True

    producer = MagicMock()
    tempfile = _make_tempfile("output.tmp")
    mock_dag.depending[producer] = {}

    assert mock_dag.is_needed_tempfile(producer, tempfile)
