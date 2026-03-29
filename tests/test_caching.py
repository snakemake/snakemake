"""
Tests for the caching module changes introduced in the feat/flag-typed PR.

Key changes tested:
- hash_file() new function in snakemake.caching.hash
- ProvenanceHashMap.get_provenance_hash and _get_provenance_hash converted from async to sync
- OutputFileCache.get_outputfiles_and_cachefiles in local.py converted from async to sync
  and changed return type from list to generator
- OutputFileCache._get_storage_objects in storage.py converted from async to sync generator
- snakemake.__init__.py version import changed to use snakemake._version module
"""

import hashlib
import inspect
import os
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Helpers to ensure the src layout is importable
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).parent.parent / "src"


def _ensure_src_in_path():
    src = str(SRC_DIR)
    if src not in sys.path:
        sys.path.insert(0, src)


_ensure_src_in_path()


# ---------------------------------------------------------------------------
# Tests for hash_file()
# ---------------------------------------------------------------------------


class TestHashFile:
    """Tests for the newly added hash_file() function in snakemake.caching.hash."""

    def _import_hash_file(self):
        from snakemake.caching.hash import hash_file

        return hash_file

    def test_hash_file_returns_hex_string(self, tmp_path):
        """hash_file should return a hex string."""
        hash_file = self._import_hash_file()
        f = tmp_path / "sample.txt"
        f.write_bytes(b"hello world")
        result = hash_file(str(f))
        assert isinstance(result, str)
        # SHA-256 hex digest is 64 characters long
        assert len(result) == 64
        # All characters must be valid hex
        int(result, 16)

    def test_hash_file_consistent_for_same_content(self, tmp_path):
        """The same content must always produce the same hash."""
        hash_file = self._import_hash_file()
        content = b"reproducible content"
        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_bytes(content)
        f2.write_bytes(content)
        assert hash_file(str(f1)) == hash_file(str(f2))

    def test_hash_file_different_content_different_hash(self, tmp_path):
        """Different file contents must produce different hashes."""
        hash_file = self._import_hash_file()
        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_bytes(b"content A")
        f2.write_bytes(b"content B")
        assert hash_file(str(f1)) != hash_file(str(f2))

    def test_hash_file_empty_file(self, tmp_path):
        """hash_file should handle empty files without errors."""
        hash_file = self._import_hash_file()
        f = tmp_path / "empty.txt"
        f.write_bytes(b"")
        result = hash_file(str(f))
        # SHA-256 of empty bytes is well-defined
        expected = hashlib.sha256(b"").hexdigest()
        assert result == expected

    def test_hash_file_matches_manual_sha256(self, tmp_path):
        """hash_file result must match a manual SHA-256 calculation."""
        hash_file = self._import_hash_file()
        content = b"test data for hashing"
        f = tmp_path / "test.dat"
        f.write_bytes(content)
        expected = hashlib.sha256(content).hexdigest()
        assert hash_file(str(f)) == expected

    def test_hash_file_large_file_block_boundary(self, tmp_path):
        """hash_file should read files in 4096-byte blocks correctly."""
        hash_file = self._import_hash_file()
        # Create a file slightly larger than one block (4096 bytes)
        content = b"x" * 8193
        f = tmp_path / "large.bin"
        f.write_bytes(content)
        expected = hashlib.sha256(content).hexdigest()
        assert hash_file(str(f)) == expected

    def test_hash_file_exactly_one_block(self, tmp_path):
        """hash_file should handle files of exactly one block size."""
        hash_file = self._import_hash_file()
        content = b"a" * 4096
        f = tmp_path / "one_block.bin"
        f.write_bytes(content)
        expected = hashlib.sha256(content).hexdigest()
        assert hash_file(str(f)) == expected

    def test_hash_file_binary_content(self, tmp_path):
        """hash_file should handle arbitrary binary content."""
        hash_file = self._import_hash_file()
        content = bytes(range(256)) * 100
        f = tmp_path / "binary.bin"
        f.write_bytes(content)
        expected = hashlib.sha256(content).hexdigest()
        assert hash_file(str(f)) == expected

    def test_hash_file_stable_across_calls(self, tmp_path):
        """Calling hash_file multiple times on the same file gives the same result."""
        hash_file = self._import_hash_file()
        f = tmp_path / "stable.txt"
        f.write_bytes(b"stable content")
        result1 = hash_file(str(f))
        result2 = hash_file(str(f))
        assert result1 == result2

    def test_hash_file_path_object(self, tmp_path):
        """hash_file should work when passed a Path object."""
        hash_file = self._import_hash_file()
        content = b"path object test"
        f = tmp_path / "pathobj.txt"
        f.write_bytes(content)
        # Passing Path directly (not str)
        result = hash_file(f)
        expected = hashlib.sha256(content).hexdigest()
        assert result == expected


# ---------------------------------------------------------------------------
# Tests for ProvenanceHashMap: verifying sync (non-async) interface
# ---------------------------------------------------------------------------


class TestProvenanceHashMapSync:
    """
    Tests that ProvenanceHashMap.get_provenance_hash and _get_provenance_hash
    are synchronous (not coroutines) after the async-to-sync conversion.
    """

    def _import_phm(self):
        from snakemake.caching.hash import ProvenanceHashMap

        return ProvenanceHashMap

    def test_get_provenance_hash_is_not_coroutine(self):
        """get_provenance_hash must not be a coroutine function (was async before PR)."""
        ProvenanceHashMap = self._import_phm()
        assert not inspect.iscoroutinefunction(ProvenanceHashMap.get_provenance_hash)

    def test_get_provenance_hash_private_is_not_coroutine(self):
        """_get_provenance_hash must not be a coroutine function (was async before PR)."""
        ProvenanceHashMap = self._import_phm()
        assert not inspect.iscoroutinefunction(ProvenanceHashMap._get_provenance_hash)

    def test_get_provenance_hash_is_regular_method(self):
        """get_provenance_hash should be a regular method, not a coroutine."""
        ProvenanceHashMap = self._import_phm()
        phm = ProvenanceHashMap()
        assert callable(phm.get_provenance_hash)
        assert not inspect.iscoroutinefunction(phm.get_provenance_hash)

    def test_provenance_hash_map_instantiation(self):
        """ProvenanceHashMap can be instantiated without errors."""
        ProvenanceHashMap = self._import_phm()
        phm = ProvenanceHashMap()
        assert phm._hashes == {}

    def test_get_provenance_hash_with_mocked_job(self, tmp_path):
        """
        get_provenance_hash should return a hex string when called synchronously
        with a properly mocked Job.
        """
        from snakemake.caching.hash import ProvenanceHashMap
        from snakemake.settings.types import DeploymentMethod

        # Build a minimal mock for a Job that has no script/shell/wrapper/notebook,
        # no input files from deps, no conda env, no container.
        job = MagicMock()
        job.rule.name = "myrule"
        job.rule.shellcmd = "echo hello"
        job.rule.script = None
        job.rule.notebook = None
        job.rule.wrapper = None
        job.is_shell = True
        job.is_script = False
        job.is_notebook = False
        job.is_wrapper = False
        job.params._allitems.return_value = []
        job.input = []
        job.dag.dependencies = {job: {}}
        job.rule.cache = MagicMock()
        job.rule.cache.omit_software = False
        job.dag.workflow.deployment_settings.deployment_method = set()
        job.conda_env = None
        job.container_img_url = None

        phm = ProvenanceHashMap()
        result = phm.get_provenance_hash(job)

        # Result must be a 64-char hex string (SHA-256)
        assert isinstance(result, str)
        assert len(result) == 64
        int(result, 16)

    def test_get_provenance_hash_is_deterministic(self):
        """
        get_provenance_hash must be deterministic for the same job description.
        """
        from snakemake.caching.hash import ProvenanceHashMap

        job = MagicMock()
        job.rule.name = "stable_rule"
        job.rule.shellcmd = "cat {input} > {output}"
        job.rule.script = None
        job.rule.notebook = None
        job.rule.wrapper = None
        job.is_shell = True
        job.is_script = False
        job.is_notebook = False
        job.is_wrapper = False
        job.params._allitems.return_value = [("key", "value")]
        job.input = []
        job.dag.dependencies = {job: {}}
        job.rule.cache = MagicMock()
        job.rule.cache.omit_software = False
        job.dag.workflow.deployment_settings.deployment_method = set()
        job.conda_env = None
        job.container_img_url = None

        phm1 = ProvenanceHashMap()
        phm2 = ProvenanceHashMap()
        assert phm1.get_provenance_hash(job) == phm2.get_provenance_hash(job)

    def test_hashes_cached_after_first_call(self):
        """_hashes dict should be populated after computing a hash."""
        from snakemake.caching.hash import ProvenanceHashMap

        job = MagicMock()
        job.rule.name = "cache_test"
        job.rule.shellcmd = "echo cached"
        job.rule.script = None
        job.rule.notebook = None
        job.rule.wrapper = None
        job.is_shell = True
        job.is_script = False
        job.is_notebook = False
        job.is_wrapper = False
        job.params._allitems.return_value = []
        job.input = []
        job.dag.dependencies = {job: {}}
        job.rule.cache = MagicMock()
        job.rule.cache.omit_software = False
        job.dag.workflow.deployment_settings.deployment_method = set()
        job.conda_env = None
        job.container_img_url = None

        phm = ProvenanceHashMap()
        assert job not in phm._hashes
        phm.get_provenance_hash(job)
        assert job in phm._hashes


# ---------------------------------------------------------------------------
# Tests for local.OutputFileCache.get_outputfiles_and_cachefiles (sync)
# ---------------------------------------------------------------------------


class TestLocalOutputFileCacheSync:
    """
    Tests for the synchronous conversion of get_outputfiles_and_cachefiles
    in snakemake.caching.local.OutputFileCache.
    """

    def _mock_jobs_module(self):
        """
        Inject a stub snakemake.jobs module so local.py can be imported
        without needing the full scheduler plugin stack.
        """
        if "snakemake.jobs" not in sys.modules:
            mock_jobs = ModuleType("snakemake.jobs")
            mock_jobs.Job = MagicMock  # placeholder
            sys.modules["snakemake.jobs"] = mock_jobs

    def _import_local_cache(self):
        self._mock_jobs_module()
        from snakemake.caching.local import OutputFileCache

        return OutputFileCache

    def test_get_outputfiles_and_cachefiles_not_coroutine(self):
        """get_outputfiles_and_cachefiles must be synchronous after the PR."""
        OutputFileCache = self._import_local_cache()
        assert not inspect.iscoroutinefunction(
            OutputFileCache.get_outputfiles_and_cachefiles
        )

    def test_get_outputfiles_and_cachefiles_returns_generator(self, tmp_path):
        """get_outputfiles_and_cachefiles should return a generator (not a list)."""
        OutputFileCache = self._import_local_cache()

        with patch.dict(os.environ, {"SNAKEMAKE_OUTPUT_CACHE": str(tmp_path)}):
            cache = OutputFileCache()
            job = MagicMock()
            job.rule.output = [MagicMock(is_multiext=False)]
            job.output = [str(tmp_path / "out.txt")]
            cache.provenance_hash_map.get_provenance_hash = MagicMock(
                return_value="abc123"
            )

            result = cache.get_outputfiles_and_cachefiles(job)
            # Must be a generator (or at least an iterator), not a list
            assert not isinstance(result, list)
            assert hasattr(result, "__iter__")
            assert hasattr(result, "__next__")

    def test_get_outputfiles_and_cachefiles_values(self, tmp_path):
        """get_outputfiles_and_cachefiles should yield (outputfile, cachefile) pairs."""
        OutputFileCache = self._import_local_cache()

        with patch.dict(os.environ, {"SNAKEMAKE_OUTPUT_CACHE": str(tmp_path)}):
            cache = OutputFileCache()
            output_path = str(tmp_path / "result.txt")
            job = MagicMock()
            job.rule.output = [MagicMock(is_multiext=False)]
            job.output = [output_path]
            fake_hash = "deadbeef" * 8  # 64 hex chars
            cache.provenance_hash_map.get_provenance_hash = MagicMock(
                return_value=fake_hash
            )

            pairs = list(cache.get_outputfiles_and_cachefiles(job))
            assert len(pairs) == 1
            out_f, cache_f = pairs[0]
            assert isinstance(out_f, Path)
            assert isinstance(cache_f, Path)
            assert out_f == Path(output_path)
            # The cache file should live under cache_location / provenance_hash + ext
            assert str(tmp_path) in str(cache_f)
            assert fake_hash in str(cache_f)
            assert cache_f.suffix == ".txt"

    def test_store_is_coroutine(self, tmp_path):
        """store() method should remain an async coroutine (it's an async method)."""
        OutputFileCache = self._import_local_cache()
        assert inspect.iscoroutinefunction(OutputFileCache.store)

    def test_fetch_is_coroutine(self, tmp_path):
        """fetch() method should remain an async coroutine."""
        OutputFileCache = self._import_local_cache()
        assert inspect.iscoroutinefunction(OutputFileCache.fetch)

    def test_exists_is_coroutine(self, tmp_path):
        """exists() method should remain an async coroutine."""
        OutputFileCache = self._import_local_cache()
        assert inspect.iscoroutinefunction(OutputFileCache.exists)


# ---------------------------------------------------------------------------
# Tests for storage.OutputFileCache._get_storage_objects (sync)
# ---------------------------------------------------------------------------


class TestStorageOutputFileCacheSync:
    """
    Tests for the synchronous conversion of _get_storage_objects
    in snakemake.caching.storage.OutputFileCache.
    """

    def _mock_deps(self):
        if "snakemake.jobs" not in sys.modules:
            mock_jobs = ModuleType("snakemake.jobs")
            mock_jobs.Job = MagicMock
            sys.modules["snakemake.jobs"] = mock_jobs

    def _import_storage_cache(self):
        self._mock_deps()
        from snakemake.caching.storage import OutputFileCache

        return OutputFileCache

    def test_get_storage_objects_not_coroutine(self):
        """_get_storage_objects must be synchronous after the PR."""
        OutputFileCache = self._import_storage_cache()
        assert not inspect.iscoroutinefunction(OutputFileCache._get_storage_objects)

    def test_get_storage_objects_not_async_generator(self):
        """_get_storage_objects must not be an async generator function."""
        OutputFileCache = self._import_storage_cache()
        assert not inspect.isasyncgenfunction(OutputFileCache._get_storage_objects)

    def test_store_is_coroutine(self):
        """store() must remain an async coroutine in storage.OutputFileCache."""
        OutputFileCache = self._import_storage_cache()
        assert inspect.iscoroutinefunction(OutputFileCache.store)

    def test_fetch_is_coroutine(self):
        """fetch() must remain an async coroutine in storage.OutputFileCache."""
        OutputFileCache = self._import_storage_cache()
        assert inspect.iscoroutinefunction(OutputFileCache.fetch)

    def test_exists_is_coroutine(self):
        """exists() must remain an async coroutine in storage.OutputFileCache."""
        OutputFileCache = self._import_storage_cache()
        assert inspect.iscoroutinefunction(OutputFileCache.exists)

    def test_get_storage_objects_yields_storage_objects(self, tmp_path):
        """_get_storage_objects should yield properly configured storage objects."""
        OutputFileCache = self._import_storage_cache()

        output_path = str(tmp_path / "output.txt")
        fake_hash = "cafebabe" * 8

        storage_object = MagicMock()
        mock_provider = MagicMock()
        mock_provider.object.return_value = storage_object

        with patch.dict(os.environ, {"SNAKEMAKE_OUTPUT_CACHE": "/cache"}):
            cache = OutputFileCache(storage_provider=mock_provider)
            cache.provenance_hash_map.get_provenance_hash = MagicMock(
                return_value=fake_hash
            )

            job = MagicMock()
            job.rule.output = [MagicMock(is_multiext=False)]
            job.output = [output_path]

            objects = list(cache._get_storage_objects(job))
            assert len(objects) == 1
            # Verify the storage provider was called with the right path
            mock_provider.object.assert_called_once()
            call_arg = mock_provider.object.call_args[0][0]
            assert fake_hash in call_arg
            # Verify set_local_path was called
            storage_object.set_local_path.assert_called_once_with(Path(output_path))

    def test_get_storage_objects_raises_on_missing_output(self, tmp_path):
        """_get_storage_objects should raise WorkflowError if output doesn't exist."""
        from snakemake.exceptions import WorkflowError

        OutputFileCache = self._import_storage_cache()

        nonexistent = str(tmp_path / "nonexistent_output.txt")
        fake_hash = "aabbccdd" * 8

        mock_provider = MagicMock()

        with patch.dict(os.environ, {"SNAKEMAKE_OUTPUT_CACHE": "/cache"}):
            cache = OutputFileCache(storage_provider=mock_provider)
            cache.provenance_hash_map.get_provenance_hash = MagicMock(
                return_value=fake_hash
            )

            job = MagicMock()
            job.rule.output = [MagicMock(is_multiext=False)]
            job.output = [nonexistent]

            with pytest.raises(WorkflowError, match="Cannot move output file"):
                list(cache._get_storage_objects(job, check_output_exists=True))

    def test_get_storage_objects_no_check_missing_output(self, tmp_path):
        """_get_storage_objects without check_output_exists should not raise for missing files."""
        OutputFileCache = self._import_storage_cache()

        nonexistent = str(tmp_path / "nonexistent_output.txt")
        fake_hash = "11223344" * 8

        storage_object = MagicMock()
        mock_provider = MagicMock()
        mock_provider.object.return_value = storage_object

        with patch.dict(os.environ, {"SNAKEMAKE_OUTPUT_CACHE": "/cache"}):
            cache = OutputFileCache(storage_provider=mock_provider)
            cache.provenance_hash_map.get_provenance_hash = MagicMock(
                return_value=fake_hash
            )

            job = MagicMock()
            job.rule.output = [MagicMock(is_multiext=False)]
            job.output = [nonexistent]

            # Should not raise
            objects = list(cache._get_storage_objects(job, check_output_exists=False))
            assert len(objects) == 1


# ---------------------------------------------------------------------------
# Tests for snakemake.__init__.py version import
# ---------------------------------------------------------------------------


class TestVersionImport:
    """Tests for the version import change in snakemake/__init__.py."""

    def test_version_is_importable(self):
        """snakemake.__version__ should be importable."""
        import snakemake

        assert hasattr(snakemake, "__version__")

    def test_version_is_string(self):
        """snakemake.__version__ must be a string."""
        import snakemake

        assert isinstance(snakemake.__version__, str)

    def test_version_not_empty(self):
        """snakemake.__version__ must not be empty."""
        import snakemake

        assert len(snakemake.__version__) > 0

    def test_version_from_version_module(self):
        """Version should be sourced from snakemake._version module."""
        from snakemake._version import version as _ver

        import snakemake

        assert snakemake.__version__ == _ver

    def test_version_module_exists(self):
        """snakemake._version module must exist and export 'version'."""
        from snakemake._version import version

        assert isinstance(version, str)