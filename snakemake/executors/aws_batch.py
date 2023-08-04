import base64
import os
import uuid
from collections import namedtuple
import time

from snakemake.common import async_lock, get_container_image
from snakemake.logging import logger
from snakemake.executors import ClusterExecutor, sleep
import heapq
import botocore
import threading
import boto3
from .aws_batch_helper import detect_aws_region, efs_id_from_access_point
from snakemake.exceptions import WorkflowError
from ..remote.S3 import S3Helper
import re

AwsBatchJob = namedtuple(
    "AwsBatchJob", "job job_id callback error_callback"
)


class AwsBatchExecutor(ClusterExecutor):
    def __init__(self, workflow, dag, cores,
                 jobname="snakejob.{name}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 restart_times=None,
                 container_image=None,
                 fsap_id=None,
                 efs_project_path=None,
                 task_queue=None,
                 workflow_role=None,
                 data_already_uploaded=False,
                 assume_shared_fs=True,
                 max_status_checks_per_second=1,
                 task_timeout=None,
                 tags=None):

        self.terminated_jobs = False
        self.max_status_checks_per_second = max_status_checks_per_second

        self.aws_region_name = detect_aws_region()
        if not self.aws_region_name:
            raise WorkflowError(
                "Failed to detect AWS region; configure AWS CLI or set environment AWS_REGION"
            )

        self.task_queue = task_queue
        self.workflow_role = workflow_role
        self.cores = cores
        self.created_job_defs = list()
        self.workflow = workflow
        self.workdir = os.path.dirname(self.workflow.persistence.path)
        self.task_timeout = task_timeout
        self.tags = tags if tags is not None else dict()

        try:
            self.aws = boto3.Session().client(  # Session() needed for thread safety
                "batch",
                region_name=self.aws_region_name,
                config=botocore.config.Config(
                    retries={"max_attempts": 5, "mode": "standard"}),
            )
        except botocore.exceptions.ClientError as exn:
            raise WorkflowError(exn)

        # this should always be True
        self.assume_shared_fs = assume_shared_fs

        if assume_shared_fs is False:
            self.s3c = S3Helper()
            self.snakemake_data = workflow.default_remote_prefix
            if data_already_uploaded is False:
                logger.debug("Uploading data to s3")
                self.upload_workflow_to_s3(workflow)
        else:
            self.fsap_id = fsap_id
            self.fs_id = efs_id_from_access_point(
                self.aws_region_name, self.fsap_id)
            self.s3_upload = False
            self.mount_path = "/mnt/efs"
            self.efs_project_path = efs_project_path

        self._describer = BatchJobDescriber()

        self.container_image = container_image or get_container_image()

        # TODO: Investigate options
        self.exec_job = (
            "snakemake {target} --snakefile %s "
            "--force --cores {cores} --keep-target-files --keep-remote "
            "--latency-wait {latency_wait} --scheduler {workflow.scheduler_type} "
            "--attempt 1 --max-inventory-time 0 "
            "{overwrite_config} {rules} --nocolor "
            "--notemp --no-hooks --nolock " % self.get_snakefile()
        )

        self._setup_job_def()

        super().__init__(workflow, dag, None,
                         jobname=jobname,
                         printreason=printreason,
                         quiet=quiet,
                         printshellcmds=printshellcmds,
                         restart_times=restart_times,
                         assume_shared_fs=assume_shared_fs,
                         max_status_checks_per_second=10)

    def upload_workflow_to_s3(self, workflow):
        s3_bucket = workflow.default_remote_prefix.split("/")[0]
        s3_subdir = re.sub(
            "^{}/".format(s3_bucket), "", workflow.default_remote_prefix
        )
        # Prepare upload the project to S3
        for root, dirs, files in os.walk(self.workdir):
            dest_path = root.replace(self.workdir, "")
            if not dest_path.startswith("/.snakemake"):
                for file in files:
                    filepath = os.path.join(root, file)
                    s3_key = "{}{}/{}".format(s3_subdir, dest_path, file)
                    logger.debug(
                        "uploading file {} to s3 location {}".format(filepath, s3_key))
                    self.s3c.upload_to_s3(s3_bucket, filepath, s3_key)

    def get_snakefile(self):
        """The snakefile must be a relative path, which should be derived
        from the self.workflow.main_snakefile.
        """
        assert os.path.exists(self.workflow.main_snakefile)
        self.workdir = os.path.dirname(self.workflow.persistence.path)
        return self.workflow.main_snakefile.replace(self.workdir, "").strip(
            os.sep
        )

    def shutdown(self):
        # perform additional steps on shutdown if necessary
        # deregister everything from AWS so the environment is clean
        logger.info("shutting down")
        if not self.terminated_jobs:
            self.terminate_active_jobs()

        for job_def in self.created_job_defs:
            self.deregister(job_def)
        super().shutdown()

    def cancel(self):
        self.terminate_active_jobs()
        self.shutdown()

    def terminate_active_jobs(self):
        for job in self.active_jobs:
            logger.debug("terminating job {}".format(job))
            try:
                self.aws.terminate_job(
                    jobId=job.job_id, reason="terminated by Snakemake")
            except botocore.exceptions.ClientError as exn:
                logger.warning(
                    "failed to terminate Batch job definition {} with error {}".format(
                        job, exn)
                )
        self.terminated_jobs = True

    def deregister(self, job_def_handle):
        try:
            logger.debug(
                "de-registering Batch job definition {}".format(job_def_handle))
            self.aws.deregister_job_definition(
                jobDefinition=job_def_handle["jobDefinitionArn"])
        except botocore.exceptions.ClientError as exn:
            # AWS expires job definitions after 6mo, so failing to delete them isn't fatal
            logger.warning(
                "failed to deregister Batch job definition {} with error {}".format(
                    job_def_handle, exn)
            )

    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):
        """Run an individual job or a job group.
        """

        super()._run(job)
        job_id = self._submit_batch_job(job)

        self.active_jobs.append(
            AwsBatchJob(job, job_id, callback, error_callback)
        )

    def _setup_job_def(self):
        # TODO: Change to use the default snakemake name
        batch_name_entropy = base64.b32encode(
            uuid.uuid4().bytes[:5]).lower().decode()
        batch_job_name = "snakemake_{}".format(batch_name_entropy)

        container_properties = {
            "executionRoleArn": self.workflow_role,
            "command": ["snakemake"],
            "privileged": True,
            "resourceRequirements": [
                {"type": "VCPU", "value": "1"},
                {"type": "MEMORY", "value": "1024"},
            ],
        }

        # TODO: Use second image and create script to execute
        container_properties["image"] = self.container_image
        container_properties["volumes"], container_properties["mountPoints"] = self._prepare_mounts()

        self.job_def = self.aws.register_job_definition(
            jobDefinitionName=batch_job_name,
            type="container",
            containerProperties=container_properties,
            tags=self.tags,
        )

        self.created_job_defs.append(self.job_def)

    def get_python_executable(self):
        return "python"

    def _generate_snakemake_command(self, job):
        exec_job = self.format_job_exec(job)

        # TODO: This needs to be a generated script that only downloads the relevant values
        command = list(filter(None, exec_job.replace("\"", "").split(" ")))
        
        return_command = ["sh", "-c"]
        snakemake_run_command = "cd {}/{} && {}".format(
            self.mount_path, self.efs_project_path, " ".join(command))
        return_command.append(snakemake_run_command)

        return return_command

    def _submit_batch_job(self, job):
        """
        Register & submit AWS batch job, leaving a cleanup callback to deregister the transient
        job definition.
        """
        # Update job name as "." is not a valid character in batch
        job_name = job.format_wildcards(self.jobname).replace(".", "_")
        command = self._generate_snakemake_command(job)

        vcpu = job.resources.get("_cores", 1)
        mem = job.resources.get("mem_mb", 1024)

        if self.task_timeout is not None:
            job = self.aws.submit_job(
                jobName=job_name,
                jobQueue=self.task_queue,
                jobDefinition="{}:{}".format(
                    self.job_def["jobDefinitionName"], self.job_def["revision"]),
                containerOverrides={
                    "command": command,
                    "resourceRequirements": [
                        {"type": "VCPU", "value": str(vcpu)},
                        {"type": "MEMORY", "value": str(mem)},
                    ],
                },
                timeout={"attemptDurationSeconds": self.task_timeout},
                tags=self.tags,
            )
        else:
            job = self.aws.submit_job(
                jobName=job_name,
                jobQueue=self.task_queue,
                jobDefinition="{}:{}".format(
                    self.job_def["jobDefinitionName"], self.job_def["revision"]),
                containerOverrides={
                    "command": command,
                    "resourceRequirements": [
                        {"type": "VCPU", "value": str(vcpu)},
                        {"type": "MEMORY", "value": str(mem)},
                    ],
                },
                tags=self.tags,
            )
        logger.debug(
            "AWS Batch job submitted with queue {}, jobId {} and tags {}".format(self.task_queue,
                                                                                 job["jobId"],
                                                                                 self.tags)
        )
        return job["jobId"]

    def _get_job_status(self, job):
        """
        Poll for Batch job success or failure & return exit code
        """
        job_id = job.job_id
        exit_code = None
        log_stream_name = None
        job_desc = self._describer.describe(self.aws, job_id, 1)
        job_status = job_desc["status"]
        if "container" in job_desc and "logStreamName" in job_desc["container"]:
            log_stream_name = job_desc["container"]["logStreamName"]
        log_details = {"status": job_status, "jobId": job_id}
        if log_stream_name:
            log_details["logStreamName"] = log_stream_name
        if job_status not in (
                "SUBMITTED",
                "PENDING",
                "RUNNABLE",
                "STARTING",
                "RUNNING",
                "SUCCEEDED",
                "FAILED",
        ):
            logger.warning(
                "unknown job status {} from AWS Batch".format(job_status))
        logger.debug("The log details is {} with status {}".format(
            log_details, job_status))
        if job_status == "SUCCEEDED":
            return 0
        elif job_status == "FAILED":
            reason = job_desc.get("container", {}).get("reason", None)
            status_reason = job_desc.get("statusReason", None)
            self.failure_info = {"jobId": job_id}
            if reason:
                self.failure_info["reason"] = reason
            if status_reason:
                self.failure_info["statusReason"] = status_reason
            if log_stream_name:
                self.failure_info["logStreamName"] = log_stream_name
            if status_reason and "Host EC2" in status_reason and "terminated" in status_reason:
                raise WorkflowError(
                    "AWS Batch job interrupted (likely spot instance termination) with error {}".format(
                        self.failure_info)
                )
            if "exitCode" not in job_desc.get("container", {}):
                raise WorkflowError(
                    "AWS Batch job failed with error {}. View log stream {} for more information".format(
                        self.failure_info["statusReason"], self.failure_info["statusReason"])
                )
            exit_code = job_desc["container"]["exitCode"]
            assert isinstance(exit_code, int) and exit_code != 0
        return exit_code

    def _prepare_mounts(self):
        """
        Prepare the "volumes" and "mountPoints" for the Batch job definition, assembling the
        in-container filesystem with the shared working directory, read-only input files, and
        command/stdout/stderr files.
        """

        # EFS mount point
        volumes = [
            {
                "name": "efs",
                "efsVolumeConfiguration": {
                    "fileSystemId": self.fs_id,
                    "transitEncryption": "ENABLED",
                },
            }
        ]
        volumes[0]["efsVolumeConfiguration"]["authorizationConfig"] = {
            "accessPointId": self.fsap_id
        }
        mount_points = [
            {"containerPath": self.mount_path, "sourceVolume": "efs"}]

        return volumes, mount_points

    async def _wait_for_jobs(self):
        # busy wait on job completion
        # This is only needed if your backend does not allow to use callbacks
        # for obtaining job status.
        while True:
            # always use self.lock to avoid race conditions
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for j in active_jobs:
                # use self.status_rate_limiter to avoid too many API calls.
                async with self.status_rate_limiter:
                    status_code = self._get_job_status(j)
                    if status_code == 0:
                        j.callback(j.job)
                    elif status_code != None:
                        j.error_callback(j.job)
                    else:
                        still_running.append(j)

            # Retrieve status of job j from your backend via j.jobid
            # Handle completion and errors, calling either j.callback(j.job)
            # or j.error_callback(j.job)
            # In case of error, add job j to still_running.
            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            await sleep()


class BatchJobDescriber:
    """
    This singleton object handles calling the AWS Batch DescribeJobs API with up to 100 job IDs
    per request, then dispensing each job description to the thread interested in it. This helps
    avoid AWS API request rate limits when we're tracking many concurrent jobs.
    """

    JOBS_PER_REQUEST = 100  # maximum jobs per DescribeJob request

    def __init__(self):
        self.lock = threading.Lock()
        self.last_request_time = 0
        self.job_queue = []
        self.jobs = {}

    def describe(self, aws, job_id, period):
        """
        Get the latest Batch job description
        """
        while True:
            with self.lock:
                if job_id not in self.jobs:
                    # register new job to be described ASAP
                    heapq.heappush(self.job_queue, (0.0, job_id))
                    self.jobs[job_id] = None
                # update as many job descriptions as possible
                self._update(aws, period)
                # return the desired job description if we have it
                desc = self.jobs[job_id]
                if desc:
                    return desc
            # otherwise wait (outside the lock) and try again
            time.sleep(period / 4)

    def unsubscribe(self, job_id):
        """
        Unsubscribe from a job_id once we'll no longer be interested in it
        """
        with self.lock:
            if job_id in self.jobs:
                del self.jobs[job_id]

    def _update(self, aws, period):
        # if enough time has passed since our last DescribeJobs request
        if time.time() - self.last_request_time >= period:
            # take the N least-recently described jobs
            job_ids = set()
            assert self.job_queue
            while self.job_queue and len(job_ids) < self.JOBS_PER_REQUEST:
                job_id = heapq.heappop(self.job_queue)[1]
                assert job_id not in job_ids
                if job_id in self.jobs:
                    job_ids.add(job_id)
            if not job_ids:
                return
            # describe them
            try:
                job_descs = aws.describe_jobs(jobs=list(job_ids))
            finally:
                # always: bump last_request_time and re-enqueue these jobs
                self.last_request_time = time.time()
                for job_id in job_ids:
                    heapq.heappush(
                        self.job_queue, (self.last_request_time, job_id))
            # update self.jobs with the new descriptions
            for job_desc in job_descs["jobs"]:
                job_ids.remove(job_desc["jobId"])
                self.jobs[job_desc["jobId"]] = job_desc
            assert not job_ids, "AWS Batch DescribeJobs didn't return all expected results"
