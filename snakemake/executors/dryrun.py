__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from snakemake_interface_executor_plugins.executors.base import AbstractExecutor
from snakemake_interface_executor_plugins.jobs import (
    JobExecutorInterface,
)
from snakemake_interface_executor_plugins import CommonSettings
from snakemake_interface_executor_plugins.executors.base import SubmittedJobInfo

from snakemake.logging import logger


common_settings = CommonSettings(
    non_local_exec=False,
    dryrun_exec=True,
    implies_no_shared_fs=False,
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
            for j in job.jobs:
                self.printcache(j)
        else:
            self.printcache(job)

    def printcache(self, job: JobExecutorInterface):
        cache_mode = self.workflow.get_cache_mode(job.rule)
        if cache_mode:
            if self.workflow.output_file_cache.exists(job, cache_mode):
                logger.info(
                    "Output file {} will be obtained from global between-workflow cache.".format(
                        job.output[0]
                    )
                )
            else:
                logger.info(
                    "Output file {} will be written to global between-workflow cache.".format(
                        job.output[0]
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
