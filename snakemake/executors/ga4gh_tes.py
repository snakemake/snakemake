__author__ = "Sven Twardziok, Alex Kanitz, Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import stat
import time
from collections import namedtuple

from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
from snakemake.executors import ClusterExecutor
from snakemake.common import Mode, get_container_image

TaskExecutionServiceJob = namedtuple(
    "TaskExecutionServiceJob", "job jobid callback error_callback"
)


class TaskExecutionServiceExecutor(ClusterExecutor):
    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob.{name}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        latency_wait=3,
        cluster_config=None,
        local_input=None,
        restart_times=None,
        assume_shared_fs=False,
        max_status_checks_per_second=0.5,
        tes_url=None,
        container_image=None,
    ):
        import tes

        self.container_image = container_image or get_container_image()
        self.container_workdir = "/tmp"
        self.max_status_checks_per_second = max_status_checks_per_second
        self.tes_url = tes_url
        self.tes_client = tes.HTTPClient(url=self.tes_url)

        logger.info("[TES] Job execution on TES: {url}".format(url=self.tes_url))

        exec_job = "\\\n".join(
            (
                "{envvars} ",
                "mkdir /tmp/conda && cd /tmp && ",
                "snakemake {target} ",
                "--snakefile {snakefile} ",
                "--verbose ",
                "--force -j{cores} ",
                "--keep-target-files ",
                "--keep-remote ",
                "--latency-wait 10 ",
                "--attempt 1 ",
                "{use_threads}",
                "{overwrite_config} {rules} ",
                "--nocolor ",
                "--notemp ",
                "--no-hooks ",
                "--nolock ",
                "--mode {} ".format(Mode.cluster),
            )
        )

        super().__init__(
            workflow,
            dag,
            None,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            cluster_config=cluster_config,
            local_input=local_input,
            restart_times=restart_times,
            exec_job=exec_job,
            assume_shared_fs=assume_shared_fs,
            max_status_checks_per_second=max_status_checks_per_second,
        )

    def write_jobscript(self, job, jobscript, **kwargs):

        use_threads = "--force-use-threads" if not job.is_group() else ""
        envvars = "\\\n".join(
            "export {}={};".format(var, os.environ[var])
            for var in self.workflow.envvars
        )

        exec_job = self.format_job(
            self.exec_job,
            job,
            _quote_all=False,
            use_threads=use_threads,
            envvars=envvars,
            **kwargs
        )
        content = self.format_job(self.jobscript, job, exec_job=exec_job, **kwargs)
        logger.debug("Jobscript:\n{}".format(content))
        with open(jobscript, "w") as f:
            print(content, file=f)
        os.chmod(jobscript, os.stat(jobscript).st_mode | stat.S_IXUSR)

    def shutdown(self):
        # perform additional steps on shutdown if necessary
        super().shutdown()

    def cancel(self):
        for job in self.active_jobs:
            try:
                self.tes_client.cancel_task(job.jobid)
                logger.info("[TES] Task canceled: {id}".format(id=job.jobid))
            except Exception:
                logger.info(
                    "[TES] Canceling task failed. This may be because the job is "
                    "already in a terminal state."
                )
        self.shutdown()

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)

        jobscript = self.get_jobscript(job)
        self.write_jobscript(job, jobscript)

        # submit job here, and obtain job ids from the backend
        try:
            task = self._get_task(job, jobscript)
            tes_id = self.tes_client.create_task(task)
            logger.info("[TES] Task submitted: {id}".format(id=tes_id))
        except Exception as e:
            raise WorkflowError(str(e))

        self.active_jobs.append(
            TaskExecutionServiceJob(job, tes_id, callback, error_callback)
        )

    def _wait_for_jobs(self):
        UNFINISHED_STATES = [
            "UNKNOWN",
            "INITIALIZING",
            "QUEUED",
            "RUNNING",
            "PAUSED",
        ]
        ERROR_STATES = [
            "EXECUTOR_ERROR",
            "SYSTEM_ERROR",
            "CANCELED",  # TODO: really call `error_callback` on this?
        ]

        while True:

            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()

            for j in active_jobs:
                with self.status_rate_limiter:  # TODO: this doesn't seem to do anything?
                    res = self.tes_client.get_task(j.jobid, view="MINIMAL")
                    logger.debug(
                        "[TES] State of task '{id}': {state}".format(
                            id=j.jobid,
                            state=res.state,
                        )
                    )
                    if res.state in UNFINISHED_STATES:
                        still_running.append(j)
                    elif res.state in ERROR_STATES:
                        logger.info("[TES] Task errored: {id}".format(id=j.jobid))
                        j.error_callback(j.job)
                    elif res.state == "COMPLETE":
                        logger.info("[TES] Task completed: {id}".format(id=j.jobid))
                        j.callback(j.job)

            with self.lock:
                self.active_jobs.extend(still_running)
            time.sleep(1 / self.max_status_checks_per_second)

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
            members_path = os.path.join(
                self.container_workdir,
                str(os.path.relpath(f)),
            )
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

    def _get_task_description(self, job):
        description = ""
        if job.is_group():
            msgs = [i.message for i in job.jobs if i.message]
            if msgs:
                description = " & ".join(msgs)
        else:
            if job.message:
                description = job.message

        return description

    def _get_task_inputs(self, job, jobscript, checkdir):
        inputs = []

        # add workflow sources to inputs
        for src in self.workflow.get_sources():
            # exclude missing, hidden, empty and build files
            if (
                not os.path.exists(src)
                or os.path.basename(src).startswith(".")
                or os.path.getsize(src) == 0
                or src.endswith(".pyc")
            ):
                continue
            inputs.append(
                self._prepare_file(
                    filename=src,
                    checkdir=checkdir,
                    pass_content=True,
                )
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
                overwrite_path=os.path.join(
                    self.container_workdir,
                    "run_snakemake.sh",
                ),
                checkdir=checkdir,
                pass_content=True,
            )
        )

        return inputs

    def _get_task_outputs(self, job, checkdir):
        outputs = []
        # add output files to outputs
        for o in job.output:
            obj = self._prepare_file(
                filename=o,
                checkdir=checkdir,
                type="Output",
            )
            if obj:
                outputs.append(obj)

        # add log files to outputs
        if job.log:
            for log in job.log:
                outputs.append(
                    self._prepare_file(
                        filename=log,
                        checkdir=checkdir,
                        type="Output",
                    )
                )

        # add benchmark files to outputs
        if hasattr(job, "benchmark") and job.benchmark:
            for benchmark in job.benchmark:
                outputs.append(
                    self._prepare_file(
                        filename=benchmark,
                        checkdir=checkdir,
                        type="Output",
                    )
                )

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

    def _get_task(self, job, jobscript):
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
        if "_cores" in job.resources:
            task["resources"]["cpu_cores"] = job.resources["_cores"]
        if "mem_mb" in job.resources:
            task["resources"]["ram_gb"] = job.resources["mem_mb"] / 1000
        if "disk_mb" in job.resources:
            task["resources"]["disk_gb"] = job.resources["disk_mb"] / 1000

        tes_task = tes.Task(**task)
        logger.debug("[TES] Built task: {task}".format(task=tes_task))
        return tes_task
