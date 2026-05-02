import sys, os, subprocess
from unittest.mock import MagicMock

from snakemake.executors import local

sys.path.insert(0, os.path.dirname(__file__))

from .common import *

from snakemake import api
from snakemake.settings import types as settings
from snakemake.spawn_jobs import SpawnedJobArgsFactory
from snakemake_interface_executor_plugins.settings import CommonSettings
import copy


def _make_common_settings(auto_deploy_default_storage_provider: bool) -> CommonSettings:
    """Create a minimal CommonSettings with auto_deploy_default_storage_provider set."""
    return CommonSettings(
        non_local_exec=True,
        implies_no_shared_fs=False,
        job_deploy_sources=False,
        auto_deploy_default_storage_provider=auto_deploy_default_storage_provider,
    )


def _make_mock_workflow(
    default_storage_provider=None,
    storage_provider_settings_keys=(),
):
    """Return a minimal mock workflow for SpawnedJobArgsFactory.precommand() testing."""
    workflow = MagicMock()
    workflow.remote_execution_settings.precommand = None
    workflow.storage_settings.default_storage_provider = default_storage_provider
    # storage_provider_settings.keys() must return the given plugin names
    workflow.storage_provider_settings.keys.return_value = list(
        storage_provider_settings_keys
    )
    return workflow


def test_precommand_auto_deploy_with_default_provider():
    """precommand includes pip install for the default storage provider."""
    workflow = _make_mock_workflow(default_storage_provider="s3")
    factory = SpawnedJobArgsFactory(workflow=workflow)
    cmd = factory.precommand(
        _make_common_settings(auto_deploy_default_storage_provider=True)
    )
    assert "snakemake-storage-plugin-s3" in cmd
    assert "pip install" in cmd


def test_precommand_auto_deploy_disabled():
    """precommand does NOT include pip install when auto_deploy is disabled."""
    workflow = _make_mock_workflow(default_storage_provider="s3")
    factory = SpawnedJobArgsFactory(workflow=workflow)
    cmd = factory.precommand(
        _make_common_settings(auto_deploy_default_storage_provider=False)
    )
    assert "pip install" not in cmd


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
