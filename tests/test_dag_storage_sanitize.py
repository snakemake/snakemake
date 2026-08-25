import asyncio
from unittest.mock import MagicMock

import pytest

# Import api first to resolve the dag <-> workflow circular import.
import snakemake.api  # noqa: F401
from snakemake.dag import DAG
from snakemake.io import AnnotatedString, IOFile


def _make_storage_iofile(path):
    """Create a real IOFile flagged as a storage output.

    The file is backed by the real local filesystem so that exists_local()
    and remove() exercise the same filesystem operations as production code,
    while the storage object itself can remain mocked.
    """
    annotated = AnnotatedString(path)
    annotated.flags = {"storage_object": MagicMock()}

    rule = MagicMock()
    rule.pathvars.apply.side_effect = lambda value: value
    rule.workflow.iocache.active = False

    return IOFile(annotated, rule=rule)


def _run(coro):
    """Run an async DAG method from a synchronous pytest test."""
    return asyncio.run(coro)


@pytest.fixture()
def mock_dag():
    """Minimal DAG state required by sanitize_local_storage_copies()."""
    dag = object.__new__(DAG)
    dag._finished = set()
    dag._needrun = set()
    dag._running = set()
    return dag


def test_sanitize_removes_local_copy_of_unfinished_nonrunning_job(tmp_path, mock_dag):
    """An unfinished, non-running job should still have its local copy removed.

    Such an output is going to be recreated by this invocation, so removing
    an existing local storage copy is the intended pre-existing behavior of
    sanitize_local_storage_copies().
    """
    local_copy = tmp_path / "output.mmi"
    local_copy.write_text("stale")

    job = MagicMock()
    job.output = [_make_storage_iofile(str(local_copy))]
    mock_dag._needrun.add(job)

    _run(mock_dag.sanitize_local_storage_copies())

    assert not local_copy.exists()


def test_sanitize_keeps_local_copy_of_running_job(tmp_path, mock_dag):
    """Do not remove storage output belonging to a currently running job.

    DAG.postprocess() can run while unrelated jobs are still executing, for
    example while resolving checkpoint dependencies. A running job may
    already have created its local storage output but not yet have completed
    Snakemake's postprocessing and been marked finished.

    Removing the local copy in that interval races with normal job completion
    handling and can make a successfully produced output appear to be missing.
    """
    local_copy = tmp_path / "output.mmi"
    local_copy.write_text("freshly produced by a running job")

    job = MagicMock()
    job.output = [_make_storage_iofile(str(local_copy))]
    mock_dag._needrun.add(job)
    mock_dag._running.add(job)

    _run(mock_dag.sanitize_local_storage_copies())

    assert local_copy.exists()


def test_sanitize_removes_local_copy_of_nonrunning_group_job(tmp_path, mock_dag):
    """A running group that doesn't contain this job must not preserve it.

    Only jobs that are members of an actually running group should be
    protected from sanitation.
    """
    local_copy = tmp_path / "output.mmi"
    local_copy.write_text("stale")

    job = MagicMock()
    job.output = [_make_storage_iofile(str(local_copy))]
    other_group = MagicMock()
    other_group.is_group.return_value = True
    other_group.__contains__.return_value = False

    mock_dag._needrun.add(job)
    mock_dag._running.add(other_group)

    _run(mock_dag.sanitize_local_storage_copies())

    assert not local_copy.exists()


def test_sanitize_keeps_local_copy_of_running_group_job(tmp_path, mock_dag):
    """Do not remove output belonging to a job inside a running group.

    The scheduler registers the GroupJob itself as running rather than each
    individual job contained in it. Therefore sanitation has to protect a job
    when its associated group is present in DAG._running.
    """
    local_copy = tmp_path / "output.mmi"
    local_copy.write_text("freshly produced by a running grouped job")

    job = MagicMock()
    job.output = [_make_storage_iofile(str(local_copy))]
    group = MagicMock()
    group.is_group.return_value = True
    group.__contains__.side_effect = lambda j: j is job

    mock_dag._needrun.add(job)
    mock_dag._running.add(group)

    _run(mock_dag.sanitize_local_storage_copies())

    assert local_copy.exists()
