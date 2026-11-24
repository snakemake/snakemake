"""Unit tests for the deployment lock and the done flag in
``snakemake.deployment``. They need no network and no conda: a fake
environment writes a directory in ``tmp_path``.
"""

import asyncio
import errno
import os
import shutil
import threading
from datetime import timedelta
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import flufl.lock._lockfile as flufl_lockfile
import pytest
from flufl.lock import Lock, TimeOutError
from snakemake_interface_common.exceptions import WorkflowError
from snakemake_interface_executor_plugins.settings import SharedFSUsage

from snakemake import deployment
from snakemake.deployment import (
    DeploymentLock,
    SoftwareDeploymentManager,
    get_deployment_done_file,
    get_deployment_lock_file,
    is_deployed,
)
from snakemake.jobs import Job


class FakeEnv:
    """The parts of a ``DeployableEnvBase`` that the deployment manager uses."""

    def __init__(self, prefix: Path, deploy_time: float = 0.0, fail: bool = False):
        self.deployment_path = prefix / "abc123"
        # Read by SoftwareDeploymentManager.prepare_envs and _collect_envs.
        self.spec = "fake env"
        self.within = None
        self.deploy_time = deploy_time
        self.fail = fail
        self.deploy_count = 0
        self.remove_count = 0
        self.active = 0
        self.max_active = 0
        self.remove_error = None
        self.during_deploy = None

    def is_deployable(self) -> bool:
        return True

    def is_cacheable(self) -> bool:
        return False

    async def managed_deploy(self) -> None:
        self.deploy_count += 1
        self.active += 1
        self.max_active = max(self.max_active, self.active)
        try:
            # Like rattler, make the directory at the start.
            self.deployment_path.mkdir()
            await asyncio.sleep(self.deploy_time)
            if self.during_deploy is not None:
                self.during_deploy()
            if self.fail:
                raise WorkflowError("Deployment of fake env failed")
            (self.deployment_path / "bin").mkdir()
        finally:
            self.active -= 1

    def managed_remove(self) -> None:
        self.remove_count += 1
        if self.remove_error is not None:
            raise self.remove_error
        shutil.rmtree(self.deployment_path)

    async def cleanup_cache(self) -> None:
        pass


def _make_manager(latency_wait: int = 0, dryrun: bool = False):
    workflow = SimpleNamespace(
        deployment_settings=SimpleNamespace(deployment_methods=set()),
        execution_settings=SimpleNamespace(latency_wait=latency_wait),
        storage_settings=SimpleNamespace(
            shared_fs_usage={SharedFSUsage.SOFTWARE_DEPLOYMENT}
        ),
        remote_exec=False,
        dryrun=dryrun,
    )
    return SoftwareDeploymentManager(workflow)


def _job(env):
    return SimpleNamespace(software_env=env)


@pytest.fixture(autouse=True)
def short_lock_times(monkeypatch):
    monkeypatch.setattr(deployment, "DEPLOYMENT_LOCK_LIFETIME", timedelta(seconds=30))
    monkeypatch.setattr(deployment, "DEPLOYMENT_LOCK_POLL", timedelta(seconds=0.2))


@pytest.fixture
def env(tmp_path):
    return FakeEnv(tmp_path)


def _complete(env):
    env.deployment_path.mkdir()
    get_deployment_done_file(env).touch()


def test_paths(env):
    assert get_deployment_lock_file(env) == env.deployment_path.parent / "abc123.lock"
    assert get_deployment_done_file(env) == env.deployment_path.parent / "abc123.done"


def test_new_env_is_deployed(env):
    asyncio.run(_make_manager().prepare_envs([_job(env)]))
    assert env.deploy_count == 1
    assert is_deployed(env)
    assert not get_deployment_lock_file(env).exists()
    # The flag is the renamed claim file. unlock() deleted the lock file.
    assert get_deployment_done_file(env).stat().st_nlink == 1
    assert sorted(env.deployment_path.parent.iterdir()) == sorted(
        [env.deployment_path, get_deployment_done_file(env)]
    )


def test_complete_env_is_not_deployed_again(env, monkeypatch):
    _complete(env)
    monkeypatch.setattr(
        deployment,
        "DeploymentLock",
        MagicMock(side_effect=AssertionError("no lock")),
    )
    asyncio.run(_make_manager().prepare_envs([_job(env)]))
    assert env.deploy_count == 0


def test_concurrent_deployments_deploy_once(env):
    env.deploy_time = 0.5
    manager = _make_manager()

    async def main():
        await asyncio.gather(manager._deploy_env(env), manager._deploy_env(env))

    asyncio.run(main())
    assert env.deploy_count == 1
    assert env.max_active == 1
    assert is_deployed(env)


def test_directory_without_flag_is_repaired(env):
    env.deployment_path.mkdir()
    (env.deployment_path / "partial").touch()
    asyncio.run(_make_manager().prepare_envs([_job(env)]))
    assert env.remove_count == 1
    assert env.deploy_count == 1
    assert is_deployed(env)
    assert not (env.deployment_path / "partial").exists()


def test_late_flag_prevents_removal(env):
    env.deployment_path.mkdir()
    timer = threading.Timer(0.3, get_deployment_done_file(env).touch)
    timer.start()
    try:
        asyncio.run(_make_manager(latency_wait=3)._deploy_env(env))
    finally:
        timer.cancel()
    assert env.remove_count == 0
    assert env.deploy_count == 0
    assert is_deployed(env)


def test_flag_without_directory_is_deployed(env):
    env.deployment_path.parent.mkdir(parents=True, exist_ok=True)
    get_deployment_done_file(env).touch()
    asyncio.run(_make_manager().prepare_envs([_job(env)]))
    assert env.deploy_count == 1
    assert env.remove_count == 0
    assert is_deployed(env)


def test_failed_deployment_writes_no_flag(env):
    env.fail = True
    with pytest.raises(WorkflowError):
        asyncio.run(_make_manager().prepare_envs([_job(env)]))
    assert env.deployment_path.exists()
    assert not get_deployment_done_file(env).exists()
    assert not get_deployment_lock_file(env).exists()

    # The next run removes the partial directory and deploys again.
    env.fail = False
    asyncio.run(_make_manager().prepare_envs([_job(env)]))
    assert env.remove_count == 1
    assert env.deploy_count == 2
    assert is_deployed(env)


def test_held_lock_logs_warning_once(env, monkeypatch):
    mock_logger = MagicMock()
    monkeypatch.setattr(deployment, "logger", mock_logger)
    env.deployment_path.parent.mkdir(parents=True, exist_ok=True)
    other = Lock(str(get_deployment_lock_file(env)), lifetime=timedelta(seconds=30))
    other.lock()
    # flufl.lock sleeps up to 2 s between tries. Hold the lock longer, so that
    # at least one try times out.
    timer = threading.Timer(4.0, other.unlock)
    timer.start()
    try:
        asyncio.run(_make_manager()._deploy_env(env))
    finally:
        timer.cancel()
    assert mock_logger.warning.call_count == 1
    assert str(get_deployment_lock_file(env)) in mock_logger.warning.call_args[0][0]
    assert env.deploy_count == 1
    assert is_deployed(env)


def test_refresher_keeps_lock_during_long_deployment(env, monkeypatch):
    # Without the refresher, a waiting task breaks the lock after 3 s.
    # flufl.lock stores the expiry time in whole seconds, so do not use a
    # lifetime shorter than this.
    monkeypatch.setattr(deployment, "DEPLOYMENT_LOCK_LIFETIME", timedelta(seconds=3))
    monkeypatch.setattr(flufl_lockfile, "CLOCK_SLOP", timedelta(0))
    env.deploy_time = 6
    manager = _make_manager()

    async def main():
        first = asyncio.create_task(manager._deploy_env(env))
        await asyncio.sleep(0.2)
        await asyncio.gather(first, manager._deploy_env(env))

    asyncio.run(main())
    assert env.deploy_count == 1
    assert env.max_active == 1
    assert is_deployed(env)


def test_lost_lock_raises_and_writes_no_flag(env):
    env.during_deploy = lambda: os.unlink(get_deployment_lock_file(env))
    with pytest.raises(WorkflowError, match="Lost deployment lock"):
        asyncio.run(_make_manager()._deploy_env(env))
    assert not get_deployment_done_file(env).exists()


def test_lock_broken_before_flag_write_writes_no_flag(env, monkeypatch):
    # Another run breaks the lock after the ownership check, just before the
    # flag is written. flufl.lock then deletes the claim file of this run, so
    # the rename of the claim file to the flag fails.
    done_file = str(get_deployment_done_file(env))
    real_replace = os.replace

    def replace(src, dst):
        if str(dst) == done_file:
            DeploymentLock(str(get_deployment_lock_file(env)))._break()
        return real_replace(src, dst)

    monkeypatch.setattr(os, "replace", replace)
    with pytest.raises(WorkflowError, match="Lost deployment lock"):
        asyncio.run(_make_manager()._deploy_env(env))
    assert not get_deployment_done_file(env).exists()
    assert not is_deployed(env)


def test_no_run_takes_lock_while_old_claim_file_exists(env, monkeypatch):
    # Run A is about to rename its claim file to the done flag. Run B breaks
    # the lock and pauses after its first delete. Run C must not get the lock
    # in this gap, because the claim file of A must already be gone.
    lock_file = str(get_deployment_lock_file(env))
    done_file = str(get_deployment_done_file(env))
    real_replace = os.replace
    real_unlink = os.unlink
    c_got_lock = []

    def unlink_then_try_c(path):
        real_unlink(path)
        if not c_got_lock:
            monkeypatch.setattr(os, "unlink", real_unlink)
            c = DeploymentLock(lock_file)
            try:
                c.lock(timedelta(seconds=0.5))
                c_got_lock.append(True)
                c.unlock()
            except TimeOutError:
                c_got_lock.append(False)

    def replace(src, dst):
        if str(dst) == done_file:
            monkeypatch.setattr(os, "unlink", unlink_then_try_c)
            DeploymentLock(lock_file)._break()
            monkeypatch.setattr(os, "unlink", real_unlink)
        return real_replace(src, dst)

    monkeypatch.setattr(os, "replace", replace)
    with pytest.raises(WorkflowError, match="Lost deployment lock"):
        asyncio.run(_make_manager()._deploy_env(env))
    assert c_got_lock == [False]
    assert not get_deployment_done_file(env).exists()


def test_no_hard_links_gives_clear_error(env, monkeypatch):
    env.deployment_path.parent.mkdir(parents=True, exist_ok=True)

    def no_link(src, dst):
        raise OSError(errno.EPERM, "Operation not permitted")

    monkeypatch.setattr(os, "link", no_link)
    with pytest.raises(WorkflowError) as excinfo:
        asyncio.run(_make_manager()._deploy_env(env))
    assert str(get_deployment_lock_file(env)) in str(excinfo.value)
    assert "hard links" in str(excinfo.value)
    assert env.deploy_count == 0


def test_dryrun_does_not_lock_remove_or_deploy(env, monkeypatch):
    env.deployment_path.mkdir()
    monkeypatch.setattr(
        deployment,
        "DeploymentLock",
        MagicMock(side_effect=AssertionError("no lock")),
    )
    asyncio.run(_make_manager(dryrun=True).prepare_envs([_job(env)]))
    assert env.deploy_count == 0
    assert env.remove_count == 0
    assert env.deployment_path.exists()
    assert not get_deployment_lock_file(env).exists()


def test_cleanup_removes_flag_when_remove_fails(env):
    _complete(env)
    env.remove_error = WorkflowError("Removal of fake env failed")
    with pytest.raises(WorkflowError):
        asyncio.run(_make_manager().cleanup_envs([_job(env)]))
    assert not get_deployment_done_file(env).exists()
    assert not get_deployment_lock_file(env).exists()


def test_cleanup_accepts_two_jobs_with_one_env(env):
    _complete(env)
    asyncio.run(_make_manager().cleanup_envs([_job(env), _job(env)]))
    assert env.remove_count == 1
    assert not env.deployment_path.exists()
    assert not get_deployment_done_file(env).exists()


def test_get_wait_for_files_returns_done_flag(env):
    job = SimpleNamespace(
        local_input=[],
        storage_input=[],
        shadow_dir=None,
        software_env=env,
        is_wrapper=False,
        dag=SimpleNamespace(
            workflow=SimpleNamespace(
                storage_settings=SimpleNamespace(
                    shared_fs_usage={SharedFSUsage.SOFTWARE_DEPLOYMENT}
                )
            )
        ),
    )
    assert Job.get_wait_for_files(job, skip_input_files=True) == [
        str(get_deployment_done_file(env))
    ]
