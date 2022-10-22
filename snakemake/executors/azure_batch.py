__author__ = "Johannes Köster, Andreas Wilm, Jake VanCampen"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
from collections import namedtuple
import datetime
import uuid
import io
import shutil
import tarfile
import tempfile
import sys

from snakemake.executors import ClusterExecutor, sleep
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.common import get_container_image, get_file_hash

# Define an Azure Batch job. Snakemake requires this namedtuple to at least
# contain the attributes: job, jobid, callback, error_callback. In AzBatch
# parlance this is actually a task and a job is the overarching container. FIXME
# should probably make self.job_id part of this as well to get rid of it
AzBatchJob = namedtuple("AzBatchJob", "job jobid task_id callback error_callback")


class AzBatchExecutor(ClusterExecutor):
    "Azure Batch Executor"

    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob.{name}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        container_image=None,
        regions=None,
        location=None,
        cache=False,
        local_input=None,
        restart_times=None,
        exec_job=None,
        max_status_checks_per_second=1,
        az_batch_account_name=None,
        az_batch_account_key=None,
        az_batch_account_location=None,
    ):


        try:
            from azure.batch import BatchServiceClient
            from azure.batch.batch_auth import SharedKeyCredentials

        except ImportError:
            raise WorkflowError(
                "The Python 3 package 'azure-batch'"
                " must be installed to use Azure Batch"
            )

        # FIXME overwrite options here or upstream? docker on, sharedfs off etc

        self.workflow = workflow
        # self.quiet = quiet
        self.workdir = os.path.dirname(self.workflow.persistence.path)
        # self._save_storage_cache = cache

        # Relative path for running on instance
        self._set_snakefile()

        # Prepare workflow sources for build package
        self._set_workflow_sources()

        # Pool ids can only contain any combination of alphanumeric characters along with dash and underscore
        ts = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        self.pool_id = "snakemakepool-{:s}".format(ts)
        self.job_id = "snakemakejob-{:s}".format(ts)

        exec_job = exec_job or (
            "snakemake {target} --snakefile %s "
            "--force -j{cores} --keep-target-files --keep-remote "
            "--attempt 1 {use_threads} "
            "{overwrite_config} {rules} --nocolor "
            "--notemp --no-hooks --nolock " % self.snakefile
        )

        self.envvars = list(self.workflow.envvars) or []

        self.container_image = container_image or get_container_image()

        # initialize az batch account information from command line args
        self.az_batch_account_name     = az_batch_account_name
        self.az_batch_account_key      = az_batch_account_key
        self.az_batch_account_location = az_batch_account_location
        self.az_batch_url              = f"{self.az_batch_account_name}.{self.az_batch_account_location}.batch.azure.com"

        # As in google_lifesciences.py:
        # Keep track of build packages to clean up shutdown, and generate
        self._build_packages = set()
        targz = self._generate_build_source_package()
        self.resource_file = self._upload_build_source_package(targz)

        # Create a Batch service client.
        credentials = SharedKeyCredentials(
            self.az_batch_account_name,
            self.az_batch_account_key,
        )
        self.batch_client = BatchServiceClient(
            credentials, batch_url=self.az_batch_url
        )

        # Create the pool that will contain the compute nodes that will execute the
        # tasks.
        logger.debug("Creating AzBatch pool %s" % self.pool_id)
        self.create_pool(self.batch_client, self.pool_id)

        # Create the job that will run the tasks.
        logger.debug("Creating AzBatch job")
        self.create_job(self.batch_client, self.job_id, self.pool_id)

        super().__init__(
            workflow,
            dag,
            None,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            restart_times=restart_times,
            exec_job=exec_job,
            assume_shared_fs=False,
            max_status_checks_per_second=10,
        )

    def shutdown(self):
        # perform additional steps on shutdown if necessary (jobs were cancelled already)
        from azure.storage.blob import BlobClient
        
        logger.debug("Deleting AzBatch job")
        self.batch_client.job.delete(self.job_id)

        logger.debug("Deleting AzBatch pool")
        self.batch_client.pool.delete(self.pool_id)

        logger.debug("Deleting uploaded workflow sources")
        blob_client = BlobClient.from_blob_url(blob_url=self.resource_file.http_url)
        blob_client.delete_blob()

        super().shutdown()

    def cancel(self):
        for task in self.batch_client.task.list(self.job_id):
            # strictly not need as job deletion also deletes task
            self.batch_client.task.terminate(self.job_id, task.id)
        self.shutdown()

    def run(self, job, callback=None, submit_callback=None, error_callback=None):

        import azure.batch._batch_service_client as batch
        import azure.batch.models as batchmodels

        super()._run(job)

        # obtain job execution command
        use_threads = "--force-use-threads" if not job.is_group() else ""

        envsettings = []
        for key in self.envvars:
            try:
                envsettings.append(
                    batchmodels.EnvironmentSetting(name=key, value=os.environ[key])
                )
            except KeyError:
                continue

        exec_job = self.format_job(
            self.exec_job, job, _quote_all=True, use_threads=use_threads
        )

        # FIXME how come we don't need source activate snakemake here? See
        # _generate_job_action() in google_lifesciences
        exec_job += self.get_default_resources_args()
        # exec_job = "/bin/sh -c 'ls && mount && pwd'"
        exec_job = "/bin/sh -c 'tar xzf {} && {}'".format(
            self.resource_file.file_path, exec_job
        )
        logger.debug("Exec job = %s" % exec_job)

        task_id = str(
            uuid.uuid1()
        )  # A string that uniquely identifies the Task within the Job.

        # useful blog https://www.muspells.net/blog/2018/11/azure-batch-task-in-containers/

        # This is the user who run the command inside the container.
        # An unprivileged one
        user = batchmodels.AutoUserSpecification(
            scope=batchmodels.AutoUserScope.task,
            # `non_admin` would be better, but conda complains
            # about non writable path  /home/_azbatchtask_1/.conda/pkgs/urls.txt
            elevation_level=batchmodels.ElevationLevel.admin,
        )

        # This is the docker image we want to run
        task_container_settings = batchmodels.TaskContainerSettings(
            image_name=self.container_image, container_run_options="--rm"
        )

        # https://docs.microsoft.com/en-us/python/api/azure-batch/azure.batch.models.taskaddparameter?view=azure-python
        # all directories recursively below the AZ_BATCH_NODE_ROOT_DIR (the root of Azure Batch directories on the node) are mapped into the container, all Task environment variables are mapped into the container, and the Task command line is executed in the container
        task = batch.models.TaskAddParameter(
            id=task_id,
            command_line=exec_job,
            container_settings=task_container_settings,
            resource_files=[self.resource_file],  # Snakefile, yml files etc.
            user_identity=batchmodels.UserIdentity(auto_user=user),
            environment_settings=envsettings,
        )

        # FIXME autorestart/retry should be handled here

        # register job as active, using your own namedtuple.
        self.batch_client.task.add(self.job_id, task)
        self.active_jobs.append(AzBatchJob(job, task_id, callback, error_callback))
        logger.debug("Added AzBatch task %s" % task_id)

    # from https://github.com/Azure-Samples/batch-python-quickstart/blob/master/src/python_quickstart_client.py
    @staticmethod
    def _read_stream_as_string(stream, encoding):
        """Read stream as string
        :param stream: input stream generator
        :param str encoding: The encoding of the file. The default is utf-8.
        :return: The file content.
        :rtype: str
        """
        output = io.BytesIO()
        try:
            for data in stream:
                output.write(data)
            if encoding is None:
                encoding = "utf-8"
            return output.getvalue().decode(encoding)
        finally:
            output.close()
        raise RuntimeError("could not write data to stream or decode bytes")

    # adopted from https://github.com/Azure-Samples/batch-python-quickstart/blob/master/src/python_quickstart_client.py
    def _get_task_output(self, job_id, task_id, stdout_or_stderr, encoding=None):
        assert stdout_or_stderr in ["stdout", "stderr"]
        fname = stdout_or_stderr + ".txt"
        stream = self.batch_client.file.get_from_task(job_id, task_id, fname)
        content = self._read_stream_as_string(stream, encoding)
        return content

    def _wait_for_jobs(self):
        # FIXME this is run early i.e self.batch_client.task.list()
        # might not return anything. so we need to keep track
        # of jobs differently. task.add doesn't return the tasks
        # just task.add
        import azure.batch.models as batchmodels

        while True:
            # always use self.lock to avoid race conditions
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()

            logger.debug("Monitoring %d active AzBatch tasks" % len(active_jobs))
            # Loop through active jobs and act on status
            for batch_job in active_jobs:
                # use self.status_rate_limiter to avoid too many API calls.
                with self.status_rate_limiter:
                    task = self.batch_client.task.get(self.job_id, batch_job.task_id)
                    # if j.task.status == batchmodels.TaskState.completed:
                    if task.state == batchmodels.TaskState.completed:
                        dt = (
                            task.execution_info.end_time
                            - task.execution_info.start_time
                        )
                        rc = task.execution_info.exit_code
                        rt = task.execution_info.retry_count
                        stderr = self._get_task_output(
                            self.job_id, batch_job.task_id, "stderr"
                        )
                        stdout = self._get_task_output(
                            self.job_id, batch_job.task_id, "stdout"
                        )
                        sys.stderr.write(
                            "FIXME task {} completed: result={} exit_code={}\n".format(
                                batch_job.task_id, task.execution_info.result, rc
                            )
                        )
                        sys.stderr.write(
                            "FIXME task {} completed: run_time={}, retry_count={}\n".format(
                                batch_job.task_id, str(dt), rt
                            )
                        )
                        sys.stderr.write(
                            "FIXME task {}: stderr='{}'\n".format(
                                batch_job.task_id, stderr
                            )
                        )
                        sys.stderr.write(
                            "FIXME task {}: stdout='{}'\n".format(
                                batch_job.task_id, stdout
                            )
                        )

                        if (
                            task.execution_info.result
                            == batchmodels.TaskExecutionResult.failure
                        ):
                            batch_job.error_callback(batch_job.job)
                        elif (
                            task.execution_info.result
                            == batchmodels.TaskExecutionResult.success
                        ):
                            batch_job.callback(batch_job.job)
                        else:
                            raise ValueError(
                                "Unknown task execution result: {}".format(
                                    task.execution_info.result
                                )
                            )

                    # The operation is still running
                    else:
                        logger.debug("Task %s not complete" % batch_job.task_id)
                        sys.stderr.write("FIXME task {}: creation_time={} state={} node_info={}\n".format(
                            batch_job.task_id, task.creation_time, task.state, task.node_info))
                        still_running.append(batch_job)

            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()

    def create_pool(self, batch_service_client, pool_id):
        """Creates a pool of compute nodes 
        """

        import azure.batch.models as batchmodels
        import azure.batch._batch_service_client as batch

        # For more information about creating pools of Linux nodes, see:
        # https://azure.microsoft.com/documentation/articles/batch-linux-nodes/

        # Specify container configuration, fetching an image
        #  https://docs.microsoft.com/en-us/azure/batch/batch-docker-container-workloads#prefetch-images-for-container-configuration
        container_config = batch.models.ContainerConfiguration(
            container_image_names=["andreaswilm/snakemaks:5.17"]
        )  # FIXME hardcoded

        image_ref = batch.models.ImageReference(
            publisher="microsoft-azure-batch",
            offer="ubuntu-server-container",
            sku="16-04-lts",
            version="latest",
        )

        # For available env vars see https://docs.microsoft.com/en-us/azure/batch/batch-compute-node-environment-variables
        # Runs in AZ_BATCH_TASK_WORKING_DIR

        # https://docs.microsoft.com/en-us/azure/batch/resource-files
        # https://docs.microsoft.com/en-us/python/api/azure-batch/azure.batch.models.resourcefile?view=azure-python
        # We could try to down the resource file (Snakefile etc) only once, i.e. per pool rather than per task using a StartTask.
        # Question is how to get it into the Docker env then.
        # workdir = "/workdir"# seems to be the default in snakemake containers?
        # rf = self.resource_file
        # start_task = batch.models.StartTask(
        #   command_line='/bin/sh -c "mkdir -p {} && tar -xvz -C {} -f {}"'.format(
        #        workdir, workdir, rf.file_path),
        #    resource_files=[rf],
        #    wait_for_success=True)
        start_task = None
        new_pool = batch.models.PoolAddParameter(
            id=pool_id,
            virtual_machine_configuration=batchmodels.VirtualMachineConfiguration(
                image_reference=image_ref,
                container_configuration=container_config,
                node_agent_sku_id="batch.node.ubuntu 16.04",
            ),  # must match image ref
            vm_size=self.az_batch_config["BATCH_POOL_VM_SIZE"],
            target_dedicated_nodes=self.az_batch_config["BATCH_POOL_NODE_COUNT"],
            target_low_priority_nodes=0,
            start_task=start_task,
        )
        batch_service_client.pool.add(new_pool)
        # FIXME max_tasks_per_node needs to be set here as well. default is 1!
        # FIXME with higher values, will task_scheduling_policy 'pack' be
        # intelligent enough to not overload node?

        # FIXME autoscaling
        # the below confiures the cluster apparently correctly, but
        # scaling doesn't happen. All values evaluate to 0
        # set target_dedicated_nodes and target_low_priority_nodes to 0
        # if using
        # https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling
        # this formula is the task-based adjustment formula taken from the Azure portal
        if False:
            formula = """// In this example, the pool size is adjusted based on the number of tasks in the queue. 
    // Note that both comments and line breaks are acceptable in formula strings.

    // Get pending tasks for the past 15 minutes.
    $samples = $ActiveTasks.GetSamplePercent(TimeInterval_Minute * 15);
    // If we have fewer than 70 percent data points, we use the last sample point, otherwise we use the maximum of last sample point and the history average.
    $tasks = $samples < 70 ? max(0, $ActiveTasks.GetSample(1)) : 
    max( $ActiveTasks.GetSample(1), avg($ActiveTasks.GetSample(TimeInterval_Minute * 15)));
    // If number of pending tasks is not 0, set targetVM to pending tasks, otherwise half of current dedicated.
    $targetVMs = $tasks > 0 ? $tasks : max(0, $TargetDedicatedNodes / 2);
    // The pool size is capped at 20, if target VM value is more than that, set it to 20. This value should be adjusted according to your use case.
    cappedPoolSize = 20;
    $TargetDedicatedNodes = max(0, min($targetVMs, cappedPoolSize));
    // Set node deallocation mode - keep nodes active only until tasks finish
    $NodeDeallocationOption = taskcompletion;"""
            response = batch_service_client.pool.enable_auto_scale(
                pool_id,
                auto_scale_formula=formula,
                auto_scale_evaluation_interval=datetime.timedelta(minutes=10),
                pool_enable_auto_scale_options=None,
                custom_headers=None,
                raw=False,
            )

    @staticmethod
    def create_job(batch_service_client, job_id, pool_id):
        """Creates a job with the specified ID, associated with the specified pool
        """
        import azure.batch._batch_service_client as batch  # https://docs.microsoft.com/en-us/azure/batch/quick-run-python uses old 6.0 API without _

        logger.debug("Creating job %s" % job_id)

        job = batch.models.JobAddParameter(
            id=job_id,
            constraints=batch.models.JobConstraints(max_task_retry_count=0),# let snakemake decide on reruns
            pool_info=batch.models.PoolInformation(pool_id=pool_id),
        )

        batch_service_client.job.add(job)

    # from google_lifesciences.py
    def _set_snakefile(self):
        """The snakefile must be a relative path, which cannot be reliably
           derived from the self.workflow.snakefile as we might have moved
           execution into a temporary directory, and the initial Snakefile
           was somewhere else on the system.
        """
        from snakemake import SNAKEFILE_CHOICES

        for snakefile in SNAKEFILE_CHOICES:
            if os.path.exists(os.path.join(self.workdir, snakefile)):
                self.snakefile = snakefile
                break

    # from google_lifesciences.py
    def _set_workflow_sources(self):
        """We only add files from the working directory that are config related
           (e.g., the Snakefile or a config.yml equivalent), or checked into git. 
        """
        self.workflow_sources = []

        for wfs in self.workflow.get_sources():
            if os.path.isdir(wfs):
                for (dirpath, dirnames, filenames) in os.walk(wfs):
                    self.workflow_sources.extend(
                        [
                            self.workflow.check_source_sizes(os.path.join(dirpath, f))
                            for f in filenames
                        ]
                    )
            else:
                self.workflow_sources.append(
                    self.workflow.check_source_sizes(os.path.abspath(wfs))
                )

    # from google_lifesciences.py
    def _generate_build_source_package(self):
        """in order for the instance to access the working directory in storage,
           we need to upload it. This file is cleaned up at the end of the run.
           We do this, and then obtain from the instance and extract.
        """
        # Workflow sources for cloud executor must all be under same workdir root
        for filename in self.workflow_sources:
            if self.workdir not in filename:
                raise WorkflowError(
                    "All source files must be present in the working directory, "
                    "{workdir} to be uploaded to a build package that respects "
                    "relative paths, but {filename} was found outside of this "
                    "directory. Please set your working directory accordingly, "
                    "and the path of your Snakefile to be relative to it.".format(
                        workdir=self.workdir, filename=filename
                    )
                )

        # We will generate a tar.gz package, renamed by hash
        tmpname = next(tempfile._get_candidate_names())
        targz = os.path.join(tempfile.gettempdir(), "snakemake-%s.tar.gz" % tmpname)
        tar = tarfile.open(targz, "w:gz")

        # Add all workflow_sources files
        for filename in self.workflow_sources:
            arcname = filename.replace(self.workdir + os.path.sep, "")
            tar.add(filename, arcname=arcname)
        logger.debug(
            "Created %s with the following contents: %s"
            % (targz, self.workflow_sources)
        )
        tar.close()

        # Rename based on hash, in case user wants to save cache
        sha256 = get_file_hash(targz)
        hash_tar = os.path.join(
            self.workflow.persistence.aux_path, "workdir-{}.tar.gz".format(sha256)
        )

        # Only copy if we don't have it yet, clean up if we do
        if not os.path.exists(hash_tar):
            shutil.move(targz, hash_tar)
        else:
            os.remove(targz)

        # We will clean these all up at shutdown
        self._build_packages.add(hash_tar)

        return hash_tar


    @staticmethod
    def _upload_build_source_package(
        targz, container_name="resources"
    ):  # FIXME hardcoded container_name
        """given a .tar.gz created for a workflow, upload it to the blob
        storage account, only if the blob doesn't already exist.
        """
        import azure.batch.models as batchmodels

        try:
            from azure.storage.blob import BlobServiceClient
            import azure.core.exceptions
        except ImportError as e:
            raise WorkflowError(
                "The Python 3 package 'azure-storage-blob' "
                "needs to be installed to use Azure Storage remote() file functionality. %s"
                % e.msg
            )

        account_url = os.getenv(
            "AZ_BLOB_ACCOUNT_URL"
        )  # same as used for storing. FIXME requires that we used a SAS
        assert (
            account_url
        ), "Need AZ_BLOB_ACCOUNT_URL for resource file and data staging, but it's not defined"
        bsc = BlobServiceClient(account_url)

        # create container if needed
        cc = bsc.get_container_client(container_name)
        try:
            cc.create_container()
        except azure.core.exceptions.ResourceExistsError:
            pass

        # upload file
        # FIXME only if it doesn't exist
        blob_name = os.path.basename(targz)
        bc = cc.get_blob_client(blob_name)
        try:
            with open(targz, "rb") as data:
                bc.upload_blob(data, blob_type="BlockBlob")
        except Exception as e:
            raise WorkflowError("Error in creating blob. %s" % str(e))

        return batchmodels.ResourceFile(http_url=bc.url, file_path=blob_name)