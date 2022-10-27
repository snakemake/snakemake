__author__ = "Johannes Köster, Andreas Wilm, Jake VanCampen"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from configparser import ConfigParser
from typing import Tuple
import os
from collections import namedtuple
from urllib.parse import urlparse 
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
from snakemake.parser import Workdir
from snakemake.remote.AzBlob import AzureStorageHelper

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
        max_status_checks_per_second=1,
        az_batch_account_url=None,
    ):
        super().__init__(
            workflow,
            dag,
            None,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            restart_times=restart_times,
            assume_shared_fs=False,
            max_status_checks_per_second=10,
        )

        try:
            from azure.batch import BatchServiceClient
            from azure.batch.batch_auth import SharedKeyCredentials

        except ImportError:
            raise WorkflowError(
                "The Python 3 package 'azure-batch'"
                " must be installed to use Azure Batch"
            )

        self.azblob_helper  = AzureStorageHelper()

        # get container from remote prefix
        self.prefix_container = str.split(workflow.default_remote_prefix, '/')[0]
        self.resources_path = self.prefix_container + "/resources"

        # setup batch configuration sets self.az_batch_config
        self.az_batch_config = self.configure_az_batch(az_batch_account_url)
        logger.debug(self.az_batch_config)

        self.workflow = workflow
        self.workdir = os.path.dirname(self.workflow.persistence.path)

        # Relative path for running on instance
        self._set_snakefile()

        # Prepare workflow sources for build package
        self._set_workflow_sources()

        # Pool ids can only contain any combination of alphanumeric characters along with dash and underscore
        ts = datetime.datetime.now().strftime("%Y-%m%dT%H-%M-%S")
        self.pool_id = "snakepool-{:s}".format(ts)
        self.job_id = "snakejob-{:s}".format(ts)

        self.envvars = list(self.workflow.envvars) or []

        self.container_image = container_image or get_container_image()

        # Package workflow sources files and upload to storage
        self._build_packages = set()
        targz = self._generate_build_source_package()
        self.resource_file = self._upload_build_source_package(targz)

        # Create a Batch service client.
        credentials = SharedKeyCredentials(
            self.az_batch_config.get("batch_account_name"),
            self.az_batch_config.get("batch_account_key"),
        )
        self.batch_client = BatchServiceClient(
            credentials, batch_url=self.az_batch_config.get("batch_account_url"),
        )

        self.create_batch_pool()
        self.create_batch_job()
       

    def shutdown(self):
        # perform additional steps on shutdown if necessary (jobs were cancelled already)

        logger.debug("Deleting AzBatch job")
        self.batch_client.job.delete(self.job_id)

        logger.debug("Deleting AzBatch pool")
        self.batch_client.pool.delete(self.pool_id)

        logger.debug("Deleting workflow sources from blob")
        self.azblob_helper.delete_from_container(self.prefix_container, self.resource_file.file_path)

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

        envsettings = []
        for key in self.envvars:
            try:
                envsettings.append(
                    batchmodels.EnvironmentSetting(name=key, value=os.environ[key])
                )
            except KeyError:
                continue

        exec_job = self.format_job_exec(job) 
        exec_job += self.get_default_resources_args()
        exec_job = f"/bin/sh -c 'tar xzf {self.resource_file.file_path} && {exec_job}'"
        logger.debug(f"Exec job = {exec_job}")

        # A string that uniquely identifies the Task within the Job.
        task_uuid = str(uuid.uuid1())
        task_id = f"{job.rule.name}-{task_uuid}" 

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
            image_name=self.container_image, container_run_options="--rm -v /:${AZ_BATCH_TASK_WORKING_DIR}"
        )

        # https://docs.microsoft.com/en-us/python/api/azure-batch/azure.batch.models.taskaddparameter?view=azure-python
        # all directories recursively below the AZ_BATCH_NODE_ROOT_DIR (the root of Azure Batch directories on the node) 
        # are mapped into the container, all Task environment variables are mapped into the container, 
        # and the Task command line is executed in the container
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
        self.active_jobs.append(AzBatchJob(job, self.job_id, task_id, callback, error_callback))
        logger.debug(f"Added AzBatch task {task_id}")

    @staticmethod
    def configure_az_batch(az_batch_account_url: str):
        """
        :param str az_batch_account_url: the azure batch account url

        _set_az_batch_config configures self.az_batch_config with the
        necessary configuraiton variables, adding default values 
        for those configuraiton variables that can be configured 
        using environment variables.
        """

        # init config dict
        config = dict()

        config["batch_account_url"] = az_batch_account_url

        
        result = urlparse(az_batch_account_url)
        batch_account_name = str.split(result.hostname,".")[0]
        config["batch_account_name"] = batch_account_name

        batch_account_key = os.getenv("AZ_BATCH_ACCOUNT_KEY")
        if batch_account_key is None:
            sys.Exit("Error: AZ_BATCH_ACCOUNT_KEY cannot be none")

        config["batch_account_key"] = batch_account_key

        # pool config defaults or env variable overrides
        config.setdefault("batch_pool_image_publisher", "microsoft-azure-batch")
        config.setdefault("batch_pool_image_offer", "ubuntu-server-container")
        config.setdefault("batch_pool_image_sku", "20-04-lts")
        config.setdefault("batch_pool_vm_container_image", "snakemake/snakemake")

        # must match batch_pool_image_sk
        config.setdefault("batch_pool_vm_node_agent_sku_id", "batch.node.ubuntu 20.04")
        config.setdefault("batch_pool_vm_size", "Standard_D2_v3")
        config.setdefault("batch_pool_node_count", 1)

        print(config)
        return config
    
    # todo cleanup config parse
    @staticmethod
    def set_config_var_if_exists(config, var):
        gotvar = os.getenv(var)
        if gotvar is not None:
            config[var] = gotvar
        return gotvar

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

    # adopted from https://github.com/Azure-Samples/batch-python-quickstart/blob/master/src/python_quickstart_client.py
    def _get_task_output(self, job_id, task_id, stdout_or_stderr, encoding=None):
        assert stdout_or_stderr in ["stdout", "stderr"]
        fname = stdout_or_stderr + ".txt"
        try:
            stream = self.batch_client.file.get_from_task(job_id, task_id, fname)
            content = self._read_stream_as_string(stream, encoding)
        except:
            content = ""
            pass


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

            logger.debug(f"Monitoring {len(active_jobs)} active AzBatch tasks")
            # Loop through active jobs and act on status
            for batch_job in active_jobs:
                task = self.batch_client.task.get(self.job_id, batch_job.task_id)

                # print(task)
                # # check compute node state to handle case where task remains active
                # # because the node is "unusable"
                # if task.node_info is not None:
                #     node = self.batch_client.compute_node.get(task.node_info.pool_id, task.node_info.node_id)
                #     print(node.state, node.errors)

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
                    logger.debug(f"Task {batch_job.task_id} state={task.state}") 
                    sys.stderr.write("FIXME task {}: creation_time={} state={} node_info={}\n".format(
                        batch_job.task_id, task.creation_time, task.state, task.node_info))
                    still_running.append(batch_job)

            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()

    def create_batch_pool(self):
        """Creates a pool of compute nodes 
        """

        import azure.batch._batch_service_client as bsc
        import azure.batch.models as batchmodels

        image_ref = bsc.models.ImageReference(
             publisher=self.az_batch_config["batch_pool_image_publisher"],
             offer=self.az_batch_config["batch_pool_image_offer"],
             sku=self.az_batch_config["batch_pool_image_sku"],
             version="latest",
        ) 

        # Specify container configuration, fetching an image
        #  https://docs.microsoft.com/en-us/azure/batch/batch-docker-container-workloads#prefetch-images-for-container-configuration
        container_config = bsc.models.ContainerConfiguration(
            container_image_names=[self.az_batch_config["batch_pool_vm_container_image"]]
            # container_image_names=['andreaswilm/snakemaks:5.17']
        ) 

        start_task = None
        new_pool = bsc.models.PoolAddParameter(
            id=self.pool_id,
            virtual_machine_configuration=batchmodels.VirtualMachineConfiguration(
                image_reference=image_ref,
                container_configuration=container_config,
                node_agent_sku_id=self.az_batch_config["batch_pool_vm_node_agent_sku_id"]
            ), 
            vm_size=self.az_batch_config["batch_pool_vm_size"],
            target_dedicated_nodes=self.az_batch_config["batch_pool_node_count"],
            target_low_priority_nodes=0,
            start_task=start_task,
        )

        # create pool if not exists
        try:
            logger.debug(f"Creating pool: {self.pool_id}")
            self.batch_client.pool.add(new_pool)
        except batchmodels.BatchErrorException as err:
            if err.error.code != "PoolExists":
                raise
            else:
                logger.debug(f"Pool {self.pool_id} exists.")


    def create_batch_job(self):
        """Creates a job with the specified ID, associated with the specified pool
        """
        import azure.batch._batch_service_client as bsc

        logger.debug(f"Creating job {self.job_id}")

        self.batch_client.job.add(
            bsc.models.JobAddParameter(
                id=self.job_id,
                constraints=bsc.models.JobConstraints(max_task_retry_count=0),
                pool_info=bsc.models.PoolInformation(pool_id=self.pool_id),
            )
        )

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

        for wfs in self.dag.get_sources():
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
                    f"All source files must be present in the working directory, "
                    "{self.workdir} to be uploaded to a build package that respects "
                    "relative paths, but {filename} was found outside of this "
                    "directory. Please set your working directory accordingly, "
                    "and the path of your Snakefile to be relative to it."
                )

        # We will generate a tar.gz package, renamed by hash
        tmpname = next(tempfile._get_candidate_names())
        targz = os.path.join(tempfile.gettempdir(), f"snakemake-{tmpname}.tar.gz")
        tar = tarfile.open(targz, "w:gz")

        # Add all workflow_sources files
        for filename in self.workflow_sources:
            arcname = filename.replace(self.workdir + os.path.sep, "")
            tar.add(filename, arcname=arcname)
        logger.debug(f"Created {targz} with the following contents: {self.workflow_sources}")
        tar.close()

        # Rename based on hash, in case user wants to save cache
        sha256 = get_file_hash(targz)
        hash_tar = os.path.join(
            self.workflow.persistence.aux_path, f"workdir-{sha256}.tar.gz"
        )

        # Only copy if we don't have it yet, clean up if we do
        if not os.path.exists(hash_tar):
            shutil.move(targz, hash_tar)
        else:
            os.remove(targz)

        # We will clean these all up at shutdown
        self._build_packages.add(hash_tar)

        return hash_tar


    def _upload_build_source_package(self, targz): 
        """given a .tar.gz created for a workflow, upload it to the blob
        storage account, only if the blob doesn't already exist.
        """

        import azure.batch.models as batchmodels
        from azure.core.exceptions import ResourceExistsError

        # create batch resources container if doesn't exist
        cc = self.azblob_helper.blob_service_client.get_container_client(self.prefix_container)
        try:
            cc.create_container()
        except ResourceExistsError:
            pass

        # upload file
        # FIXME only if it doesn't exist
        blob_name = "resources/" + os.path.basename(targz)
        bc = cc.get_blob_client(blob_name)
        try:
            with open(targz, "rb") as data:
                bc.upload_blob(data, blob_type="BlockBlob")
        except Exception as e:
            raise WorkflowError(f"Error in creating blob. {str(e)}")

        return batchmodels.ResourceFile(http_url=bc.url, file_path=blob_name)