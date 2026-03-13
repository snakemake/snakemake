__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from typing import cast

from snakemake_interface_executor_plugins.executors.base import (
    AbstractExecutor,
    SubmittedJobInfo,
)
from snakemake_interface_executor_plugins.jobs import (
    GroupJobExecutorInterface,
    JobExecutorInterface,
    SingleJobExecutorInterface,
)
from snakemake_interface_executor_plugins.settings import CommonSettings

from snakemake.logging import logger

common_settings = CommonSettings(
    non_local_exec=False,
    dryrun_exec=True,
    implies_no_shared_fs=False,
    job_deploy_sources=False,
    pass_envvar_declarations_to_cmd=False,
    auto_deploy_default_storage_provider=False,
)


class Executor(AbstractExecutor):
    def run_job(
        self,
        job: JobExecutorInterface,
    ):
        job_info = SubmittedJobInfo(job=job)
        self.report_job_submission(job_info)
        self.report_job_success(job_info)

    def get_exec_mode(self):
        raise NotImplementedError()

    def printjob(self, job: JobExecutorInterface):
        super().printjob(job)
        if job.is_group():
            group_job = cast(GroupJobExecutorInterface, job)
            for j in group_job.jobs:
                self.printcache(j)
            return
        self.printcache(cast(SingleJobExecutorInterface, job))

    def printcache(self, job: SingleJobExecutorInterface):
        if (
            self.workflow.workflow_settings.cache is not None
            and job.rule.cache
            and job.rule.cache.output
        ):
            logger.info(
                "Output file {} will be {} global between-workflow cache.".format(
                    job.output[0],
                    (
                        "obtained from"
                        if self.workflow.async_run(job.rule.cache.exists(job))
                        else "written to"
                    ),
                )
            )

    def cancel(self):
        # nothing to do
        pass

    def shutdown(self):
        # nothing to do
        pass

    def handle_job_success(self, job: JobExecutorInterface):
        # nothing to do
        pass

    def handle_job_error(self, job: JobExecutorInterface):
        # nothing to do
        pass

    @property
    def cores(self):
        return self.workflow.resource_settings.cores
