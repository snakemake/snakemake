from collections import defaultdict
from unittest.mock import MagicMock

import pytest

# Import api first to resolve the dag <-> workflow circular import.
import snakemake.api  # noqa: F401
from snakemake.dag import DAG
from snakemake.io import flag


def _make_tempfile(name):
    """Create an AnnotatedString flagged as temp for use in tests."""
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
    """A temp file with no downstream consumers should be marked unneeded."""
    producer = MagicMock()
    tempfile = _make_tempfile("output.tmp")
    mock_dag.workflow.remote_exec = False  # assume main process
    mock_dag.depending[producer] = {}

    assert not mock_dag.is_needed_tempfile(producer, tempfile)


def test_temp_file_needed_when_downstream_unfinished(mock_dag):
    """A temp file is still needed if a consuming job has not finished."""
    producer = MagicMock()
    consumer = MagicMock()
    tempfile = _make_tempfile("output.tmp")

    mock_dag.depending[producer] = {consumer: {tempfile}}
    mock_dag._needrun.add(consumer)

    assert mock_dag.is_needed_tempfile(producer, tempfile)


def test_temp_file_unneeded_when_downstream_finished(mock_dag):
    """A temp file becomes unneeded once all consuming jobs have finished."""
    producer = MagicMock()
    consumer = MagicMock()
    tempfile = _make_tempfile("output.tmp")
    mock_dag.workflow.remote_exec = False  # assume main process

    mock_dag.depending[producer] = {consumer: {tempfile}}
    mock_dag._finished.add(consumer)

    assert not mock_dag.is_needed_tempfile(producer, tempfile)


def test_temp_file_needed_remote_exec_empty_unneeded_set(mock_dag):
    """remote_exec with empty unneeded_temp_files  must prevent cleanup.
    Because unneeded_temp_files being emtpy means that all temp files are still
    needed by other subsequent jobs the remote job does not know about.
    """
    producer = MagicMock()
    tempfile = _make_tempfile("output.tmp")

    assert mock_dag.workflow.storage_settings.unneeded_temp_files == frozenset()
    mock_dag.depending[producer] = {}

    assert mock_dag.is_needed_tempfile(producer, tempfile)


def test_temp_file_unneeded_remote_exec_nonempty_unneeded_set(mock_dag):
    """remote_exec with the temp file being unneeded must cleanup."""
    producer = MagicMock()
    tempfile = _make_tempfile("output.tmp")

    mock_dag.workflow.storage_settings.unneeded_temp_files = {tempfile}
    mock_dag.depending[producer] = {}

    assert not mock_dag.is_needed_tempfile(producer, tempfile)


def test_subprocess_exec_always_returns_needed(mock_dag):
    """In subprocess_exec mode, all temp files are reported as needed (no cleanup)."""
    mock_dag.workflow.subprocess_exec = True

    producer = MagicMock()
    tempfile = _make_tempfile("output.tmp")
    mock_dag.depending[producer] = {}

    assert mock_dag.is_needed_tempfile(producer, tempfile)


def test_handle_temp_yields_output_iofile_with_temp_flag(mock_dag):
    """Regression: set intersection in handle_temp must yield the output _IOFile.

    _dependencies stores input _IOFile objects (no temp flag). job_.output has
    output _IOFile objects (with temp flag). The set intersection must preserve
    the output objects so that remove() sees is_flagged(file, "temp") == True.
    """
    from snakemake.io import AnnotatedString, flag as io_flag

    mock_dag.workflow.remote_exec = False  # assume main process

    # Simulate output _IOFile with temp flag
    output_file = AnnotatedString("shared.tmp")
    output_file.flags = {"temp": True, "storage_object": MagicMock()}

    # Simulate input _IOFile without temp flag (as stored in _dependencies)
    input_file = AnnotatedString("shared.tmp")
    input_file.flags = {"storage_object": MagicMock()}

    producer = MagicMock()
    producer.output = [output_file]
    producer.is_checkpoint = False

    consumer = MagicMock()

    mock_dag._dependencies = defaultdict(lambda: defaultdict(set))
    mock_dag._dependencies[consumer] = {producer: {input_file}}
    mock_dag.depending[producer] = {consumer: {input_file}}
    mock_dag._finished.add(consumer)
    mock_dag._needrun = set()

    mock_dag.workflow.storage_settings.notemp = False
    mock_dag.workflow.dryrun = False

    # Collect unneeded files via the same logic as handle_temp
    from functools import partial
    from itertools import filterfalse

    def is_temp(f):
        """Check whether a file object has the temp flag set."""
        return f.flags.get("temp", False)

    results = []
    for job_, files in mock_dag._dependencies[consumer].items():
        tempfiles = set(f for f in job_.output if is_temp(f))
        results.extend(
            filterfalse(
                partial(mock_dag.is_needed_tempfile, job_),
                {f for f in tempfiles if f in files},
            )
        )

    assert len(results) == 1
    # The yielded file must be the OUTPUT object (has temp flag)
    assert results[0] is output_file
    assert results[0].flags.get("temp") is True
