import sys, os, subprocess

from snakemake.executors import local

sys.path.insert(0, os.path.dirname(__file__))

from .common import *

from snakemake import api
from snakemake.settings import types as settings
import copy


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
