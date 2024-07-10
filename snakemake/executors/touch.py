__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import time

from snakemake_interface_executor_plugins.executors.real import RealExecutor
from snakemake_interface_executor_plugins.dag import DAGExecutorInterface
from snakemake_interface_executor_plugins.workflow import WorkflowExecutorInterface
from snakemake_interface_executor_plugins.logging import LoggerExecutorInterface
from snakemake_interface_executor_plugins.jobs import (
    JobExecutorInterface,
)
from snakemake_interface_executor_plugins.executors.base import SubmittedJobInfo
from snakemake_interface_executor_plugins.settings import CommonSettings
from snakemake.common import async_run

from snakemake.exceptions import print_exception


common_settings = CommonSettings(
    non_local_exec=False,
    implies_no_shared_fs=False,
    job_deploy_sources=False,
    touch_exec=True,
    pass_envvar_declarations_to_cmd=False,
    auto_deploy_default_storage_provider=False,
)


class Executor(RealExecutor):
    def run_job(
        self,
        job: JobExecutorInterface,
    ):
        job_info = SubmittedJobInfo(job=job)
        try:
            time.sleep(0.1)

            if job.output:

                async def touch():
                    for f in job.output:
                        if f.is_storage and await f.exists_in_storage():
                            await f.touch_storage_and_local()
                        elif await f.exists_local():
                            f.touch()

                async_run(touch())

            self.report_job_submission(job_info)
            self.report_job_success(job_info)
        except OSError as ex:
            print_exception(ex, self.workflow.linemaps)
            self.report_job_error(job_info)

    def get_exec_mode(self):
        raise NotImplementedError()

    def handle_job_success(self, job: JobExecutorInterface):
        super().handle_job_success(job)

    def cancel(self):
        # nothing to do
        pass

    def shutdown(self):
        # nothing to do
        pass

    def get_python_executable(self):
        raise NotImplementedError()

    @property
    def cores(self):
        return self.workflow.resource_settings.cores
