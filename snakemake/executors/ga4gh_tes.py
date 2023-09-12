__author__ = "Sven Twardziok, Alex Kanitz, Valentin Schneider-Lunitz, Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import asyncio
import math
import os
from collections import namedtuple

from snakemake_interface_executor_plugins.dag import DAGExecutorInterface
from snakemake_interface_executor_plugins.jobs import ExecutorJobInterface
from snakemake_interface_executor_plugins.workflow import WorkflowExecutorInterface
from snakemake_interface_executor_plugins.utils import sleep
from snakemake_interface_executor_plugins.executors.remote import RemoteExecutor
from snakemake_interface_executor_plugins.persistence import StatsExecutorInterface
from snakemake_interface_executor_plugins.logging import LoggerExecutorInterface

from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
from snakemake.common import get_container_image, async_lock

TaskExecutionServiceJob = namedtuple(
    "TaskExecutionServiceJob", "job jobid callback error_callback"
)


class TaskExecutionServiceExecutor(RemoteExecutor):
    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
        jobname="snakejob.{name}.{jobid}.sh",
        max_status_checks_per_second=0.5,
        tes_url=None,
        container_image=None,
    ):
        super().__init__(
            workflow,
            dag,
            stats,
            logger,
            None,
            jobname=jobname,
            max_status_checks_per_second=max_status_checks_per_second,
        )
        try:
            import tes
        except ImportError:
            raise WorkflowError(
                "Unable to import Python package tes. TES backend requires py-tes to be installed. Please install py-tes, e.g. via Conda or Pip."
            )

        self.container_image = container_image or get_container_image()
        logger.info(f"Using {self.container_image} for TES jobs.")
        self.container_workdir = "/tmp"
        self.max_status_checks_per_second = max_status_checks_per_second
        self.tes_url = tes_url
        self.tes_client = tes.HTTPClient(
            url=self.tes_url,
            token=os.environ.get("TES_TOKEN"),
            user=os.environ.get("FUNNEL_SERVER_USER"),
            password=os.environ.get("FUNNEL_SERVER_PASSWORD"),
        )
        logger.info(f"[TES] Job execution on TES: {self.tes_url}")

    def get_job_exec_prefix(self, job: ExecutorJobInterface):
        return "mkdir /tmp/conda && cd /tmp"

    def shutdown(self):
        # perform additional steps on shutdown if necessary
        super().shutdown()

    def cancel(self):
        for job in self.active_jobs:
            try:
                self.tes_client.cancel_task(job.jobid)
                logger.info(f"[TES] Task canceled: {job.jobid}")
            except Exception:
                logger.info(
                    "[TES] Canceling task failed. This may be because the job is "
                    "already in a terminal state."
                )
        self.shutdown()

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
        super()._run(job)

        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        # submit job here, and obtain job ids from the backend
        try:
            task = self._get_task(job, jobscript)
            tes_id = self.tes_client.create_task(task)
            logger.info(f"[TES] Task submitted: {tes_id}")
        except Exception as e:
            raise WorkflowError(str(e))

        self.active_jobs.append(
            TaskExecutionServiceJob(job, tes_id, callback, error_callback)
        )

    async def _wait_for_jobs(self):
        UNFINISHED_STATES = ["UNKNOWN", "INITIALIZING", "QUEUED", "RUNNING", "PAUSED"]
        ERROR_STATES = [
            "EXECUTOR_ERROR",
            "SYSTEM_ERROR",
            "CANCELED",  # TODO: really call `error_callback` on this?
        ]

        while True:
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()

            for j in active_jobs:
                async with self.status_rate_limiter:  # TODO: this doesn't seem to do anything?
                    res = self.tes_client.get_task(j.jobid, view="MINIMAL")
                    logger.debug(
                        "[TES] State of task '{id}': {state}".format(
                            id=j.jobid, state=res.state
                        )
                    )
                    if res.state in UNFINISHED_STATES:
                        still_running.append(j)
                    elif res.state in ERROR_STATES:
                        logger.info(f"[TES] Task errored: {j.jobid}")
                        j.error_callback(j.job)
                    elif res.state == "COMPLETE":
                        logger.info(f"[TES] Task completed: {j.jobid}")
                        j.callback(j.job)

            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            await asyncio.sleep(1 / self.max_status_checks_per_second)

    def _check_file_in_dir(self, checkdir, f):
        if checkdir:
            checkdir = checkdir.rstrip("/")
            if not f.startswith(checkdir):
                direrrmsg = (
                    "All files including Snakefile, "
                    + "conda env files, rule script files, output files "
                    + "must be in the same working directory: {} vs {}"
                )
                raise WorkflowError(direrrmsg.format(checkdir, f))

    def _get_members_path(self, overwrite_path, f):
        if overwrite_path:
            members_path = overwrite_path
        else:
            members_path = os.path.join(self.container_workdir, str(os.path.relpath(f)))
        return members_path

    def _prepare_file(
        self,
        filename,
        overwrite_path=None,
        checkdir=None,
        pass_content=False,
        type="Input",
    ):
        import tes

        # TODO: handle FTP files
        max_file_size = 131072
        if type not in ["Input", "Output"]:
            raise ValueError("Value for 'model' has to be either 'Input' or 'Output'.")

        members = {}

        # Handle remote files
        if hasattr(filename, "is_remote") and filename.is_remote:
            return None

        # Handle local files
        else:
            f = os.path.abspath(filename)

            self._check_file_in_dir(checkdir, f)

            members["path"] = self._get_members_path(overwrite_path, f)

            members["url"] = "file://" + f
            if pass_content:
                source_file_size = os.path.getsize(f)
                if source_file_size > max_file_size:
                    logger.warning(
                        "Will not pass file '{f}' by content, as it exceeds the "
                        "minimum supported file size of {max_file_size} bytes "
                        "defined in the TES specification. Will try to upload "
                        "file instead.".format(f=f, max_file_size=max_file_size)
                    )
                else:
                    with open(f) as stream:
                        members["content"] = stream.read()
                    members["url"] = None

        model = getattr(tes.models, type)
        logger.warning(members)
        return model(**members)

    def _get_task_description(self, job: ExecutorJobInterface):
        description = ""
        if job.is_group():
            msgs = [i.message for i in job.jobs if i.message]
            if msgs:
                description = " & ".join(msgs)
        else:
            if job.message:
                description = job.message

        return description

    def _get_task_inputs(self, job: ExecutorJobInterface, jobscript, checkdir):
        inputs = []

        # add workflow sources to inputs
        for src in self.dag.get_sources():
            # exclude missing, hidden, empty and build files
            if (
                not os.path.exists(src)
                or os.path.basename(src).startswith(".")
                or os.path.getsize(src) == 0
                or src.endswith(".pyc")
            ):
                continue
            inputs.append(
                self._prepare_file(filename=src, checkdir=checkdir, pass_content=True)
            )

        # add input files to inputs
        for i in job.input:
            obj = self._prepare_file(filename=i, checkdir=checkdir)
            if obj:
                inputs.append(obj)

        # add jobscript to inputs
        inputs.append(
            self._prepare_file(
                filename=jobscript,
                overwrite_path=os.path.join(self.container_workdir, "run_snakemake.sh"),
                checkdir=checkdir,
                pass_content=True,
            )
        )

        return inputs

    def _append_task_outputs(self, outputs, files, checkdir):
        for file in files:
            obj = self._prepare_file(filename=file, checkdir=checkdir, type="Output")
            if obj:
                outputs.append(obj)
        return outputs

    def _get_task_outputs(self, job: ExecutorJobInterface, checkdir):
        outputs = []
        # add output files to outputs
        outputs = self._append_task_outputs(outputs, job.output, checkdir)

        # add log files to outputs
        if job.log:
            outputs = self._append_task_outputs(outputs, job.log, checkdir)

        # add benchmark files to outputs
        if hasattr(job, "benchmark") and job.benchmark:
            outputs = self._append_task_outputs(outputs, job.benchmark, checkdir)

        return outputs

    def _get_task_executors(self):
        import tes

        executors = []
        executors.append(
            tes.models.Executor(
                image=self.container_image,
                command=[  # TODO: info about what is executed is opaque
                    "/bin/bash",
                    os.path.join(self.container_workdir, "run_snakemake.sh"),
                ],
                workdir=self.container_workdir,
            )
        )
        return executors

    def _get_task(self, job: ExecutorJobInterface, jobscript):
        import tes

        checkdir, _ = os.path.split(self.snakefile)

        task = {}
        task["name"] = job.format_wildcards(self.jobname)
        task["description"] = self._get_task_description(job)
        task["inputs"] = self._get_task_inputs(job, jobscript, checkdir)
        task["outputs"] = self._get_task_outputs(job, checkdir)
        task["executors"] = self._get_task_executors()
        task["resources"] = tes.models.Resources()

        # define resources
        if job.resources.get("_cores") is not None:
            task["resources"].cpu_cores = job.resources["_cores"]
        if job.resources.get("mem_mb") is not None:
            task["resources"].ram_gb = math.ceil(job.resources["mem_mb"] / 1000)
        if job.resources.get("disk_mb") is not None:
            task["resources"].disk_gb = math.ceil(job.resources["disk_mb"] / 1000)

        tes_task = tes.Task(**task)
        logger.debug(f"[TES] Built task: {tes_task}")
        return tes_task
