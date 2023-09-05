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
    ExecutorJobInterface,
)
from snakemake_interface_executor_plugins import CommonSettings

from snakemake.exceptions import print_exception


common_settings = CommonSettings(
    non_local_exec=False,
)


class Executor(RealExecutor):
    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        logger: LoggerExecutorInterface,
    ):
        super().__init__(
            workflow,
            dag,
            logger,
            executor_settings=None,
            pass_envvar_declarations_to_cmd=False,
        )

    def run_job(
        self,
        job: ExecutorJobInterface,
    ):
        try:
            # Touching of output files will be done by handle_job_success
            time.sleep(0.1)
            self.report_job_success(job)
        except OSError as ex:
            print_exception(ex, self.workflow.linemaps)
            self.report_job_error(job)

    def get_exec_mode(self):
        raise NotImplementedError()

    def handle_job_success(self, job: ExecutorJobInterface):
        super().handle_job_success(job, ignore_missing_output=True)

    def cancel(self):
        pass

    def shutdown(self):
        pass

    def get_python_executable(self):
        raise NotImplementedError()

    @property
    def cores(self):
        return self.workflow.resource_settings.cores
