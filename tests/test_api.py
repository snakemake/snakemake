import sys, os, subprocess
from unittest.mock import MagicMock, patch

from snakemake.executors import local

sys.path.insert(0, os.path.dirname(__file__))

from .common import *

from snakemake import api
from snakemake.api import ApiError
from snakemake.settings import types as settings
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
from snakemake_interface_executor_plugins.settings import CommonSettings, SharedFSUsage
import copy
import pytest


def test_deploy_sources(s3_storage):
    s3_prefix, s3_settings = s3_storage

    with api.SnakemakeApi(
        settings.OutputSettings(
            verbose=True,
            show_failed_logs=True,
        ),
    ) as snakemake_api:
        workflow_api = snakemake_api.workflow(
            storage_settings=settings.StorageSettings(
                default_storage_prefix=s3_prefix,
                default_storage_provider="s3",
                shared_fs_usage=frozenset(),
            ),
            resource_settings=settings.ResourceSettings(
                cores=1,
            ),
            storage_provider_settings=s3_settings,
            snakefile=Path(dpath("test_deploy_sources/Snakefile")),
        )
        dag_api = workflow_api.dag()

        workflow = dag_api.workflow_api._workflow
        # add dummy remote execution settings as we do not actually execute here
        # (in reality they are present)
        workflow.remote_execution_settings = settings.RemoteExecutionSettings()
        workflow._prepare_dag(
            forceall=False,
            ignore_incomplete=False,
            lock_warn_only=False,
        )
        workflow._build_dag()
        workflow.upload_sources()

        cmd = workflow.spawned_job_args_factory.precommand(local.common_settings)
        assert cmd

        origdir = os.getcwd()
        env = copy.copy(os.environ)
        env.update(workflow.spawned_job_args_factory.envvars())
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                subprocess.run(cmd, shell=True, check=True, env=env)
            finally:
                os.chdir(origdir)

# ---------------------------------------------------------------------------
# Helpers for TestDryrunExecutorValidation
# ---------------------------------------------------------------------------

# A minimal CommonSettings that represents a *remote* executor which handles
# its own file staging (can_transfer_local_files=True).
_FAKE_REMOTE_COMMON_SETTINGS = CommonSettings(
    non_local_exec=True,  # local_exec == False
    implies_no_shared_fs=False,
    job_deploy_sources=False,
    can_transfer_local_files=True,  # lets shared_fs_usage=frozenset() pass validation
)


def _make_fake_remote_plugin() -> MagicMock:
    """Return a mock executor plugin whose CommonSettings represents a remote
    executor that handles its own file staging (can_transfer_local_files=True)."""
    plugin = MagicMock()
    plugin.common_settings = _FAKE_REMOTE_COMMON_SETTINGS
    return plugin


@pytest.fixture()
def patch_registry_with_fake_remote():
    """Patch ExecutorPluginRegistry.get_plugin so that the sentinel name
    ``"fake_remote"`` resolves to the mock plugin, while all other names
    (e.g. ``"dryrun"``, ``"touch"``, ``"local"``) fall through to the real
    registry.

    This allows tests to inject a controlled CommonSettings without relying on
    any particular third-party executor plugin being installed.
    """
    fake_plugin = _make_fake_remote_plugin()

    # Capture the real function object before patch replaces it on the class.
    _real_get_plugin = ExecutorPluginRegistry.get_plugin
    _registry = ExecutorPluginRegistry()

    def _get_plugin(name: str):
        if name == "fake_remote":
            return fake_plugin
        # Call the real implementation directly, bypassing the patch.
        return _real_get_plugin(_registry, name)

    with patch.object(ExecutorPluginRegistry, "get_plugin", side_effect=_get_plugin):
        yield


class TestDryrunExecutorValidation:
    """Regression tests for issue #3973.

    When ``--dry-run`` or ``--touch`` is used, the ``execution_executor``
    parameter is set to ``"dryrun"`` or ``"touch"`` to control *how* jobs are
    executed (i.e. not at all, or by just touching outputs). However,
    workflow validation must still run against the *intended* executor
    (``executor``) so that its CommonSettings are checked — not those of the
    dryrun/touch stand-in.

    The key point is that the validation executor and the execution executor
    can have different CommonSettings. For example, an executor that handles
    its own file staging (``can_transfer_local_files=True``) should pass
    ``shared_fs_usage=frozenset()`` validation regardless of whether
    ``--dry-run`` is active. The dryrun and touch executors themselves have
    ``can_transfer_local_files=False``, so if they were mistakenly used for
    validation the check would fail.

    The ``patch_registry_with_fake_remote`` fixture injects a controlled
    CommonSettings under the sentinel name ``"fake_remote"``, so the tests do
    not depend on any specific third-party executor plugin being installed.
    """

    @staticmethod
    def _snakefile():
        return Path(dpath("test_dryrun_executor_validation/Snakefile"))

    def test_dryrun_with_remote_executor_no_shared_fs(
        self, patch_registry_with_fake_remote
    ):
        """Dry-run with an executor whose CommonSettings has
        can_transfer_local_files=True should accept shared_fs_usage=frozenset()
        without requiring a default storage provider.

        Before the fix this raised ApiError because validation ran against the
        dryrun executor's CommonSettings (can_transfer_local_files=False)
        instead of the intended executor's CommonSettings
        (can_transfer_local_files=True).
        """
        with api.SnakemakeApi(
            settings.OutputSettings(verbose=True, show_failed_logs=True),
        ) as snakemake_api:
            workflow_api = snakemake_api.workflow(
                resource_settings=settings.ResourceSettings(cores=1, nodes=3),
                storage_settings=settings.StorageSettings(
                    shared_fs_usage=frozenset(),
                ),
                snakefile=self._snakefile(),
            )
            dag_api = workflow_api.dag()
            # Should NOT raise: validation uses fake_remote's CommonSettings
            # (can_transfer_local_files=True), not the dryrun stand-in's.
            dag_api.execute_workflow(
                executor="fake_remote",
                execution_executor="dryrun",
            )

    def test_touch_with_remote_executor_no_shared_fs(
        self, patch_registry_with_fake_remote
    ):
        """Touch with an executor whose CommonSettings has
        can_transfer_local_files=True should accept shared_fs_usage=frozenset(),
        for the same reason as the dry-run case above.

        workflow.execute is mocked to avoid downstream touch-execution side
        effects; only the validation block in execute_workflow() is under test.
        """
        with api.SnakemakeApi(
            settings.OutputSettings(verbose=True, show_failed_logs=True),
        ) as snakemake_api:
            workflow_api = snakemake_api.workflow(
                resource_settings=settings.ResourceSettings(cores=1, nodes=3),
                storage_settings=settings.StorageSettings(
                    shared_fs_usage=frozenset(),
                ),
                snakefile=self._snakefile(),
            )
            dag_api = workflow_api.dag()
            # Mock workflow.execute to avoid downstream touch-execution side
            # effects; we only care that execute_workflow() validation passes.
            with patch.object(
                dag_api.workflow_api._workflow, "execute"
            ) as mock_execute:
                dag_api.execute_workflow(
                    executor="fake_remote",
                    execution_executor="touch",
                )
                mock_execute.assert_called_once()

    def test_dryrun_without_execution_executor_rejects_no_shared_fs(self):
        """When the dryrun executor is passed as executor= without an
        execution_executor override, its own CommonSettings govern validation.
        Because the dryrun executor has can_transfer_local_files=False,
        shared_fs_usage=frozenset() must still be rejected.

        This also verifies the greedy-scheduler fix: because execution_executor
        defaults to executor when not supplied, the scheduler optimisation fires
        correctly for direct executor="dryrun" API calls too.
        """
        with api.SnakemakeApi(
            settings.OutputSettings(verbose=True, show_failed_logs=True),
        ) as snakemake_api:
            workflow_api = snakemake_api.workflow(
                resource_settings=settings.ResourceSettings(cores=1),
                storage_settings=settings.StorageSettings(
                    shared_fs_usage=frozenset(),
                ),
                snakefile=self._snakefile(),
            )
            dag_api = workflow_api.dag()
            with pytest.raises(ApiError, match="default storage provider"):
                dag_api.execute_workflow(executor="dryrun")
