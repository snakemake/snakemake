__authors__ = ["Yun Jiang"]
__copyright__ = "Copyright 2026, Yun Jiang"
__license__ = "MIT"

import tempfile
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, os.path.dirname(__file__))

from snakemake.persistence import RECORD_FORMAT_VERSION, ParamsChange

from snakemake.persistence.file import FilePersistence


"""Tests for snakemake.persistence — backup/restore, code caching,
locking, and record format."""


def _make_persistence(tmp_path, nolock=False):
    mock_dag = MagicMock()
    mock_dag.workflow.dag_settings.max_checksum_file_size = 100000
    mock_dag.output_files = []
    return FilePersistence(
        nolock=nolock,
        dag=mock_dag,
        path=tmp_path / ".snakemake",
    )


class TestBackupRestore:
    """Tests for output backup and restore."""

    def test_backup_and_restore(self, tmp_path):
        p = _make_persistence(tmp_path)
        output_file = tmp_path / "result.txt"
        output_file.write_text("hello world")

        p.backup_output(output_file)
        assert output_file.exists()

        output_file.unlink()
        success = p.restore_output(output_file)
        assert success
        assert output_file.read_text() == "hello world"

    def test_restore_nonexistent_backup(self, tmp_path):
        p = _make_persistence(tmp_path)
        output_file = tmp_path / "missing.txt"
        assert not p.restore_output(output_file)

    def test_cleanup_backup(self, tmp_path):
        p = _make_persistence(tmp_path)
        output_file = tmp_path / "result.txt"
        output_file.write_text("data")
        p.backup_output(output_file)
        p.cleanup_backup(output_file)
        assert not p._get_backup_path(output_file).exists()

    def test_backup_directory(self, tmp_path):
        p = _make_persistence(tmp_path)
        output_dir = tmp_path / "result_dir"
        output_dir.mkdir()
        (output_dir / "file.txt").write_text("content")

        p.backup_output(output_dir)
        assert output_dir.exists()

        import shutil

        shutil.rmtree(output_dir)
        assert not output_dir.exists()

        success = p.restore_output(output_dir)
        assert success
        assert (output_dir / "file.txt").read_text() == "content"


class TestCodeCaching:
    """Tests for _code() method — rule source code hashing."""

    def test_code_for_shell_rule(self, tmp_path):
        p = _make_persistence(tmp_path)
        mock_rule = MagicMock()
        mock_rule.shellcmd = "echo hello"
        mock_rule.run_func_src = None
        assert p._code(mock_rule) == "echo hello"

    def test_code_for_run_rule(self, tmp_path):
        p = _make_persistence(tmp_path)
        mock_rule = MagicMock()
        mock_rule.shellcmd = None
        mock_rule.run_func_src = "def run(input, output): pass"
        assert p._code(mock_rule) == "def run(input, output): pass"

    def test_code_for_script_rule_returns_none(self, tmp_path):
        p = _make_persistence(tmp_path)
        mock_rule = MagicMock()
        mock_rule.shellcmd = None
        mock_rule.run_func_src = None
        assert p._code(mock_rule) is None

    def test_code_is_cached(self, tmp_path):
        p = _make_persistence(tmp_path)
        mock_rule = MagicMock()
        mock_rule.shellcmd = "bwa mem"
        mock_rule.run_func_src = None
        # Call twice — lru_cache should return same object
        assert p._code(mock_rule) is p._code(mock_rule)


class TestInputLogParams:
    """Tests for _input(), _log(), _params() helpers."""

    def test_log_sorted(self, tmp_path):
        p = _make_persistence(tmp_path)
        mock_job = MagicMock()
        mock_job.log = ["b.log", "a.log", "c.log"]
        assert p._log(mock_job) == ["a.log", "b.log", "c.log"]

    def test_output_sorted(self, tmp_path):
        p = _make_persistence(tmp_path)
        mock_job = MagicMock()
        mock_job.output = ["z.out", "a.out"]
        assert p._output(mock_job) == ["a.out", "z.out"]


class TestLocking:
    """Tests for lock/unlock mechanism."""

    def test_lock_and_unlock(self, tmp_path):
        mock_dag = MagicMock()
        mock_dag.workflow.dag_settings.max_checksum_file_size = 100000
        mock_dag.output_files = []
        mock_dag.jobs = []
        mock_dag.needrun_jobs = lambda: []

        p = FilePersistence(
            dag=mock_dag,
            path=tmp_path / ".snakemake",
        )

        with p.lock():
            # Lock files are created during lock, but with empty jobs
            # p.locked checks for pre-existing lock files from other processes
            pass
        # After unlock, lockfile should be cleaned up
        lockdir = tmp_path / ".snakemake" / "locks"
        lock_files = list(lockdir.glob("*.lock"))
        assert len(lock_files) == 0

    def test_nolock_is_noop(self, tmp_path):
        p = _make_persistence(tmp_path, nolock=True)
        with p.lock():
            pass  # should not raise

    def test_warn_only_lock(self, tmp_path):
        mock_dag = MagicMock()
        mock_dag.workflow.dag_settings.max_checksum_file_size = 100000
        mock_dag.output_files = []
        mock_dag.jobs = []
        mock_dag.needrun_jobs = lambda: []

        p = FilePersistence(
            warn_only=True,
            dag=mock_dag,
            path=tmp_path / ".snakemake",
        )
        with p.lock():
            pass  # should not raise even if locked


class TestB64Encoding:
    """Tests for _b64id base64 encoding."""

    def test_b64id_deterministic(self, tmp_path):
        p = _make_persistence(tmp_path)
        assert p._b64id("test.txt") == p._b64id("test.txt")

    def test_b64id_different_for_different_paths(self, tmp_path):
        p = _make_persistence(tmp_path)
        assert p._b64id("a.txt") != p._b64id("b.txt")


class TestParamsChange:
    """Tests for ParamsChange dataclass."""

    def test_empty_params_change_is_falsy(self):
        pc = ParamsChange()
        assert not pc

    def test_nonempty_params_change_is_truthy(self):
        pc = ParamsChange(only_old={"x"}, only_new=set(), files={"out.txt"})
        assert pc

    def test_or_combines_changes(self):
        a = ParamsChange(only_old={"x"}, only_new=set(), files={"f1.txt"})
        b = ParamsChange(only_old=set(), only_new={"y"}, files={"f2.txt"})
        combined = a | b
        assert combined.only_old == {"x"}
        assert combined.only_new == {"y"}
        assert combined.files == {"f1.txt", "f2.txt"}

    def test_or_with_empty(self):
        empty = ParamsChange()
        real = ParamsChange(only_old={"x"}, only_new=set(), files={"f.txt"})
        assert (empty | real) == real
        assert (real | empty) == real

    def test_iter_yields_files(self):
        pc = ParamsChange(only_old={"x"}, only_new=set(), files={"f1.txt", "f2.txt"})
        assert set(pc) == {"f1.txt", "f2.txt"}

    def test_str_empty(self):
        pc = ParamsChange()
        assert "No params change" in str(pc)

    def test_str_nonempty(self):
        pc = ParamsChange(
            only_old={"old_param"}, only_new={"new_param"}, files={"f.txt"}
        )
        s = str(pc)
        assert "old_param" in s
        assert "new_param" in s


class TestRecordFormatVersion:
    """Ensure RECORD_FORMAT_VERSION is current."""

    def test_version_is_non_decreasing(self):
        assert isinstance(RECORD_FORMAT_VERSION, int)
        assert RECORD_FORMAT_VERSION >= 6


class TestDirectoryStructure:
    """Tests that Persistence creates expected directory structure."""

    def test_creates_dot_snakemake_dirs(self, tmp_path):
        p = _make_persistence(tmp_path)
        base = tmp_path / ".snakemake"
        # These are auto-created in __init__
        assert (base / "metadata").is_dir()
        assert (base / "incomplete").is_dir()
        assert (base / "shadow").is_dir()
        assert (base / "auxiliary").is_dir()
        assert (base / "iocache").is_dir()
        assert (base / "locks").is_dir()

    def test_path_property(self, tmp_path):
        p = _make_persistence(tmp_path)
        assert p.path == tmp_path / ".snakemake"

    def test_aux_path_property(self, tmp_path):
        p = _make_persistence(tmp_path)
        assert p.aux_path == tmp_path / ".snakemake" / "auxiliary"


class TestDeactivateCache:
    """Tests for cache deactivation."""

    def test_deactivate_cache(self, tmp_path):
        p = _make_persistence(tmp_path)
        # Should not raise
        p.deactivate_cache()
        # After deactivation, _read_record should be uncached
        assert p._read_record == p._read_record_uncached
        assert p._incomplete_cache is False


class TestDropIOCache:
    def test_drop_iocache_missing_file_is_noop(self):
        """Regression for race when parallel workers concurrently drop the iocache.

        The previous check-then-remove pattern raised FileNotFoundError if a
        sibling worker deleted the file between os.path.exists() and os.remove().
        drop_iocache() must tolerate the file already being gone.
        """
        with tempfile.TemporaryDirectory() as tmpdirname:
            persistence = FilePersistence(
                dag=MagicMock(), path=Path(tmpdirname) / ".snakemake"
            )
            # iocache file does not exist; drop_iocache must not raise.
            assert not Path(persistence._iocache_filename).exists()
            persistence.drop_iocache()

    def test_drop_iocache_removes_existing_file(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            persistence = FilePersistence(
                dag=MagicMock(), path=Path(tmpdirname) / ".snakemake"
            )
            iocache_file = Path(persistence._iocache_filename)
            iocache_file.write_bytes(b"placeholder")
            assert iocache_file.exists()

            persistence.drop_iocache()
            assert not iocache_file.exists()
