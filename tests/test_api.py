import sys, os, subprocess

from snakemake.executors import local

sys.path.insert(0, os.path.dirname(__file__))

from .common import *

from snakemake import api
from snakemake.api import ApiError
from snakemake.settings import types as settings
from snakemake_interface_executor_plugins.settings import SharedFSUsage
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


class TestDryrunExecutorValidation:
    """Regression tests for issue #3973.

    When --dry-run or --touch is used, the execution_executor override should
    be "dryrun" or "touch", but validation must still run against the *intended*
    executor so that its CommonSettings (e.g. can_transfer_local_files) are
    checked, not those of the dryrun/touch stand-in.
    """

    @staticmethod
    def _snakefile():
        return Path(dpath("test_dryrun_executor_validation/Snakefile"))

    def test_dryrun_with_remote_executor_no_shared_fs(self):
        """Dry-run with a remote executor that can transfer local files should
        accept --shared-fs-usage none (empty frozenset) without requiring a
        default storage provider.

        Before the fix, this raised ApiError because validation ran against the
        dryrun executor plugin (can_transfer_local_files=False) instead of the
        intended executor (htcondor, can_transfer_local_files=True).
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
            # Should NOT raise — htcondor has can_transfer_local_files=True
            dag_api.execute_workflow(
                executor="htcondor",
                execution_executor="dryrun",
            )

    def test_touch_with_remote_executor_no_shared_fs(self):
        """Touch with a remote executor that can transfer local files should
        accept --shared-fs-usage none, just like dry-run.

        We verify that the validation in execute_workflow() passes by mocking
        workflow.execute so we don't run into downstream touch execution logic
        (which has its own separate constraints around source deployment).
        """
        from unittest.mock import patch

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
            # Mock workflow.execute to avoid downstream execution side effects;
            # we only care that execute_workflow() validation does not raise.
            with patch.object(
                dag_api.workflow_api._workflow, "execute"
            ) as mock_execute:
                dag_api.execute_workflow(
                    executor="htcondor",
                    execution_executor="touch",
                )
                mock_execute.assert_called_once()

    def test_dryrun_without_execution_executor_rejects_no_shared_fs(self):
        """When execution_executor is NOT used (old behaviour), passing the
        dryrun executor directly with --shared-fs-usage none should still be
        rejected, since dryrun has can_transfer_local_files=False."""
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
