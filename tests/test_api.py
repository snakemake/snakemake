import sys

sys.path.insert(0, os.path.dirname(__file__))

from .common import *

from snakemake import api
from snakemake import settings


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
            ),
            storage_provider_settings=s3_settings,
            snakefile=dpath("test01/Snakefile"),
        )

        workflow = workflow_api._get_workflow(check_envvars=False)
        workflow._prepare_dag()
        workflow._build_dag()
        workflow.upload_sources()

        cmd = workflow.spawned_job_args_factory.precommand()

        origdir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                subprocess.run(cmd, shell=True, check=True)
            finally:
                os.chdir(origdir)
