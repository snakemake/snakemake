import asyncio
import os
from unittest.mock import AsyncMock, MagicMock

import pytest

# Import api first to resolve the dag <-> workflow circular import.
import snakemake.api  # noqa: F401
from snakemake.dag import DAG
from snakemake.io import AnnotatedString, IOCache, IOFile, _IOFile


def _make_storage_iofile(path, workflow):
    """Create a real IOFile flagged as a storage output.

    The file is backed by the real local filesystem so that exists_local()
    and remove() exercise the same filesystem operations as production code,
    while the storage object itself can remain mocked.
    """
    annotated = AnnotatedString(path)
    annotated.flags = {"storage_object": MagicMock()}

    rule = MagicMock()
    rule.pathvars.apply.side_effect = lambda value: value
    rule.workflow = workflow

    return IOFile(annotated, rule=rule)


def _run(coro):
    """Run an async DAG method from a synchronous pytest test."""
    return asyncio.run(coro)


@pytest.fixture()
def mock_dag():
    """Minimal DAG state required by sanitize_local_storage_copies()."""
    dag = object.__new__(DAG)
    dag.workflow = MagicMock()
    dag.workflow.iocache = IOCache(max_wait_time=0)
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
    job.output = [_make_storage_iofile(str(local_copy), mock_dag.workflow)]
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
    job.output = [_make_storage_iofile(str(local_copy), mock_dag.workflow)]
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
    job.output = [_make_storage_iofile(str(local_copy), mock_dag.workflow)]
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
    job.output = [_make_storage_iofile(str(local_copy), mock_dag.workflow)]
    group = MagicMock()
    group.is_group.return_value = True
    group.__contains__.side_effect = lambda j: j is job

    mock_dag._needrun.add(job)
    mock_dag._running.add(group)

    _run(mock_dag.sanitize_local_storage_copies())

    assert local_copy.exists()


@pytest.mark.parametrize("cache_active", [True, False])
def test_sanitize_invalidates_cached_local_existence(tmp_path, mock_dag, cache_active):
    local_copy = tmp_path / "output.txt"
    local_copy.write_text("stale download")
    output = _make_storage_iofile(str(local_copy), mock_dag.workflow)
    job = MagicMock(output=[output])
    mock_dag._needrun.add(job)

    if not cache_active:
        mock_dag.workflow.iocache.deactivate()
    assert _run(output.exists_local())
    if cache_active:
        assert mock_dag.workflow.iocache.exists_local[output]
    _run(mock_dag.sanitize_local_storage_copies())

    assert not local_copy.exists()
    assert not _run(output.exists_local())
    if cache_active:
        assert mock_dag.workflow.iocache.exists_local[output] is False
    else:
        assert output not in mock_dag.workflow.iocache.exists_local


def test_sanitize_preserves_write_protected_output(tmp_path, mock_dag):
    local_copy = tmp_path / "protected.txt"
    local_copy.write_text("completed protected output")
    output = _make_storage_iofile(str(local_copy), mock_dag.workflow)
    local_copy.chmod(0o444)
    job = MagicMock(output=[output])
    mock_dag._needrun.add(job)
    try:
        if os.access(local_copy, os.W_OK):
            pytest.skip("Current user can write to read-only files (e.g. root)")
        assert _run(output.protected())
        _run(mock_dag.sanitize_local_storage_copies())
        assert local_copy.read_text() == "completed protected output"
        assert _run(output.exists_local())
    finally:
        local_copy.chmod(0o644)


def test_sanitize_removes_writable_download_even_for_protected_rule(tmp_path, mock_dag):
    local_copy = tmp_path / "retained-download.txt"
    local_copy.write_text("stale download")
    output = _make_storage_iofile(str(local_copy), mock_dag.workflow)
    output.file.flags["protected"] = True
    job = MagicMock(output=[output])
    mock_dag._needrun.add(job)

    assert not _run(output.protected())
    _run(mock_dag.sanitize_local_storage_copies())

    assert not local_copy.exists()
    assert not _run(output.exists_local())


def test_failed_sanitize_does_not_cache_output_as_absent(
    tmp_path, mock_dag, monkeypatch
):
    local_copy = tmp_path / "output.txt"
    local_copy.write_text("download")
    output = _make_storage_iofile(str(local_copy), mock_dag.workflow)
    job = MagicMock(output=[output])
    mock_dag._needrun.add(job)
    remove = AsyncMock(side_effect=PermissionError("cannot remove local copy"))
    monkeypatch.setattr(_IOFile, "remove", remove)

    assert _run(output.exists_local())
    with pytest.raises(ExceptionGroup) as error:
        _run(mock_dag.sanitize_local_storage_copies())

    assert isinstance(error.value.exceptions[0], PermissionError)
    remove.assert_awaited_once_with(remove_non_empty_dir=True, only_local=True)
    assert local_copy.exists()
    assert _run(output.exists_local())
