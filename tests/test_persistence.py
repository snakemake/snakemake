import asyncio
import tempfile
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from snakemake.persistence import MetadataRecord
from snakemake.persistence.db import DbPersistence
from snakemake.persistence.file import FilePersistence


class TestCleanupContainers:
    @patch("snakemake.deployment.singularity.Image", autospec=True, create=True)
    @patch("snakemake.dag.DAG", autospec=True, create=True)
    def test_one_unrequired_container_gets_removed(self, mock_dag, mock_img):
        snakecode = """rule all:
    input:
        "foo.txt"

rule foo:
    output:
        "foo.txt"
    container:
        "docker://quay.io/mbhall88/rasusa:0.7.0"
    shell:
        "rasusa --help &> {output}"
"""
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirpath = Path(tmpdirname)
            singularity_dir = tmpdirpath / ".snakemake/singularity"
            singularity_dir.mkdir(parents=True)
            snakefile = tmpdirpath / "Snakefile"
            snakefile.write_text(snakecode)

            unrequired_img_path = (
                singularity_dir / "32077ccc4d05977ef8b94ee6f74073fd.simg"
            )
            unrequired_img_path.touch()
            assert unrequired_img_path.exists()

            required_img_path = (
                singularity_dir / "0581af3d3099c1cc8cf0088c8efe1439.simg"
            )
            required_img_path.touch()
            assert required_img_path.exists()
            mock_img.path = str(required_img_path)
            container_imgs = {"docker://quay.io/mbhall88/rasusa:0.7.0": mock_img}
            mock_dag.container_imgs = container_imgs
            persistence = FilePersistence(dag=mock_dag)
            persistence.container_img_path = str(singularity_dir)

            persistence.cleanup_containers()

            assert required_img_path.exists()
            assert not unrequired_img_path.exists()


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


def _make_dag_mock():
    dag = MagicMock()
    dag.workflow.execution_settings.latency_wait = 5
    return dag


def _patch_get_key(persistence):
    persistence._get_key = lambda f: f


class _MockOutput(str):
    """str subclass that also supports ``await f.mtime()`` and ``await f.exists()``."""

    def __new__(cls, name, mtime=None):
        instance = super().__new__(cls, name)
        instance._mtime = mtime
        return instance

    async def mtime(self):
        result = MagicMock()
        result.local_or_storage.return_value = self._mtime or time.time()
        return result

    async def exists(self):
        return True


def _make_mock_job(outputs):
    job = MagicMock()
    job.output = outputs
    return job


def _make_finished_job(outputs):
    job = MagicMock()
    job.output = [_MockOutput(name) for name in outputs]
    job.input = []
    job.rule.name = "test_rule"
    job.rule.shellcmd = "echo test"
    job.rule.run_func_src = None
    job.shellcmd = "echo test"
    job.container_img_url = None
    job.conda_env = None
    job.log = []
    job.non_derived_params = []
    job.env_modules = None
    return job


@pytest.fixture(params=["file", "db"])
def persistence(request, tmp_path):
    if request.param == "file":
        p = FilePersistence(dag=MagicMock(), path=tmp_path / ".snakemake")
    else:
        p = DbPersistence(dag=_make_dag_mock(), path=tmp_path / ".snakemake")
    _patch_get_key(p)
    p._code = lambda rule: "test_code"
    p._input = lambda job: []
    p._log = lambda job: []
    p._params = lambda job: []
    p._conda_env = lambda job: None
    p._software_stack_hash = lambda job: "hash"
    return p


def tick() -> float:
    """Sleep briefly to guarantee a clock tick, then return current time."""
    time.sleep(0.002)
    return time.time()


class TestStarttimeRecording:
    def test_incomplete_marker_after_started(self, persistence):
        before = tick()
        persistence.started(_make_mock_job(["out.txt"]))
        after = tick()

        st = persistence._get_recorded_starttime("out.txt")
        assert st is not None
        assert before <= st <= after

    def test_rerun_replaces_stale_starttime(self, persistence):
        rec = MetadataRecord(rule="foo", starttime=1000.0, endtime=2000.0)
        persistence._write_record("out.txt", rec)
        assert persistence._read_record("out.txt").starttime == 1000.0

        before = tick()
        persistence.started(_make_mock_job(["out.txt"]))
        after = tick()

        st = persistence._get_recorded_starttime("out.txt")
        assert st is not None
        assert st > 1000.0
        assert before <= st <= after

    def test_normal_completion(self, persistence):
        before = tick()
        persistence.started(_make_mock_job(["out.txt"]))
        asyncio.run(persistence.finished(_make_finished_job(["out.txt"])))
        after = tick()

        persistence._clear_cache()
        rec = persistence._read_record("out.txt")
        assert rec is not None
        assert rec.starttime is not None
        assert before <= rec.starttime <= after
        assert rec.endtime >= rec.starttime
        assert rec.incomplete is False

    def test_finished_twice_preserves_starttime(self, persistence):
        persistence.started(_make_mock_job(["out.txt"]))
        asyncio.run(persistence.finished(_make_finished_job(["out.txt"])))
        persistence._clear_cache()
        rec1 = persistence._read_record("out.txt")
        first_starttime = rec1.starttime

        asyncio.run(persistence.finished(_make_finished_job(["out.txt"])))
        persistence._clear_cache()
        rec2 = persistence._read_record("out.txt")
        assert rec2.starttime == first_starttime

    def test_finished_without_started_fallback(self, persistence):
        # call finished directly without calling started first, should have fallback time for starttime
        asyncio.run(persistence.finished(_make_finished_job(["out.txt"])))
        persistence._clear_cache()

        rec = persistence._read_record("out.txt")
        assert rec is not None
        assert rec.starttime is not None
