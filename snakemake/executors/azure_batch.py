__author__ = "Johannes Köster, Andreas Wilm, Jake VanCampen"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import datetime
import io
import os
import re
import shutil
import shlex
import sys
import tarfile
import tempfile
import uuid
from collections import namedtuple
from pprint import pformat
from typing import Optional
from urllib.parse import urlparse

from snakemake_interface_executor_plugins.executors.remote import RemoteExecutor
from snakemake_interface_executor_plugins.persistence import StatsExecutorInterface
from snakemake_interface_executor_plugins.logging import LoggerExecutorInterface
from snakemake_interface_executor_plugins.dag import DAGExecutorInterface
from snakemake_interface_executor_plugins.jobs import ExecutorJobInterface
from snakemake_interface_executor_plugins.workflow import WorkflowExecutorInterface
from snakemake_interface_executor_plugins.utils import sleep

from snakemake.exceptions import WorkflowError
import msrest.authentication as msa

from snakemake.common import async_lock, bytesto, get_container_image, get_file_hash
from snakemake.exceptions import WorkflowError
from snakemake.executors import sleep
from snakemake.logging import logger

AzBatchJob = namedtuple("AzBatchJob", "job jobid task_id callback error_callback")


def check_source_size(filename, warning_size_gb=0.2):
    """A helper function to check the filesize, and return the file
    to the calling function Additionally, given that we encourage these
    packages to be small, we set a warning at 200MB (0.2GB).
    """
    gb = bytesto(os.stat(filename).st_size, "g")
    if gb > warning_size_gb:
        logger.warning(
            f"File {filename} (size {gb} GB) is greater than the {warning_size_gb} GB "
            f"suggested size. Consider uploading larger files to storage first."
        )
    return filename


class AzBatchConfig:
    def __init__(self, batch_account_url: str):
        # configure defaults
        self.batch_account_url = batch_account_url

        # parse batch account name
        result = urlparse(self.batch_account_url)
        self.batch_account_name = str.split(result.hostname, ".")[0]

        self.batch_account_key = self.set_or_default("AZ_BATCH_ACCOUNT_KEY", None)

        # optional subnet config
        self.batch_pool_subnet_id = self.set_or_default("BATCH_POOL_SUBNET_ID", None)

        # managed identity resource id configuration
        self.managed_identity_resource_id = self.set_or_default(
            "BATCH_MANAGED_IDENTITY_RESOURCE_ID", None
        )

        # parse subscription and resource id
        if self.managed_identity_resource_id is not None:
            self.subscription_id = self.managed_identity_resource_id.split("/")[2]
            self.resource_group = self.managed_identity_resource_id.split("/")[4]

        self.managed_identity_client_id = self.set_or_default(
            "BATCH_MANAGED_IDENTITY_CLIENT_ID", None
        )

        if self.batch_pool_subnet_id is not None:
            if (
                self.managed_identity_client_id is None
                or self.managed_identity_resource_id is None
            ):
                sys.exit(
                    "Error: BATCH_MANAGED_IDENTITY_RESOURCE_ID, BATCH_MANAGED_IDENTITY_CLIENT_ID must be set when deploying batch nodes into a private subnet!"
                )

            # parse account details necessary for batch client authentication steps
            if self.batch_pool_subnet_id.split("/")[2] != self.subscription_id:
                raise WorkflowError(
                    "Error: managed identity must be in the same subscription as the batch pool subnet."
                )

            if self.batch_pool_subnet_id.split("/")[4] != self.resource_group:
                raise WorkflowError(
                    "Error: managed identity must be in the same resource group as the batch pool subnet."
                )

        # sas url to a batch node start task bash script
        self.batch_node_start_task_sasurl = os.getenv("BATCH_NODE_START_TASK_SAS_URL")

        # options configured with env vars or default
        self.batch_pool_image_publisher = self.set_or_default(
            "BATCH_POOL_IMAGE_PUBLISHER", "microsoft-azure-batch"
        )
        self.batch_pool_image_offer = self.set_or_default(
            "BATCH_POOL_IMAGE_OFFER", "ubuntu-server-container"
        )
        self.batch_pool_image_sku = self.set_or_default(
            "BATCH_POOL_IMAGE_SKU", "20-04-lts"
        )
        self.batch_pool_vm_container_image = self.set_or_default(
            "BATCH_POOL_VM_CONTAINER_IMAGE", "ubuntu"
        )
        self.batch_pool_vm_node_agent_sku_id = self.set_or_default(
            "BATCH_POOL_VM_NODE_AGENT_SKU_ID", "batch.node.ubuntu 20.04"
        )
        self.batch_pool_vm_size = self.set_or_default(
            "BATCH_POOL_VM_SIZE", "Standard_D2_v3"
        )

        # dedicated pool node count
        self.batch_pool_node_count = self.set_or_default("BATCH_POOL_NODE_COUNT", 1)

        # default tasks per node
        # see https://learn.microsoft.com/en-us/azure/batch/batch-parallel-node-tasks
        self.batch_tasks_per_node = self.set_or_default("BATCH_TASKS_PER_NODE", 1)

        # possible values "spread" or "pack"
        # see https://learn.microsoft.com/en-us/azure/batch/batch-parallel-node-tasks
        self.batch_node_fill_type = self.set_or_default(
            "BATCH_NODE_FILL_TYPE", "spread"
        )

        # enables simplified batch node communication if set
        # see: https://learn.microsoft.com/en-us/azure/batch/simplified-compute-node-communication
        self.batch_node_communication_mode = self.set_or_default(
            "BATCH_NODE_COMMUNICATION_SIMPLIFIED", None
        )

        self.resource_file_prefix = self.set_or_default(
            "BATCH_POOL_RESOURCE_FILE_PREFIX", "resource-files"
        )

        self.container_registry_url = self.set_or_default(
            "BATCH_CONTAINER_REGISTRY_URL", None
        )

        self.container_registry_user = self.set_or_default(
            "BATCH_CONTAINER_REGISTRY_USER", None
        )

        self.container_registry_pass = self.set_or_default(
            "BATCH_CONTAINER_REGISTRY_PASS", None
        )

    @staticmethod
    def set_or_default(evar: str, default: Optional[str]):
        gotvar = os.getenv(evar)
        if gotvar is not None:
            return gotvar
        else:
            return default


# the usage of this credential helper is required to authenitcate batch with managed identity credentials
# because not all Azure SDKs support the azure.identity credentials yet, and batch is one of them.
# ref1: https://gist.github.com/lmazuel/cc683d82ea1d7b40208de7c9fc8de59d
# ref2: https://gist.github.com/lmazuel/cc683d82ea1d7b40208de7c9fc8de59d
class AzureIdentityCredentialAdapter(msa.BasicTokenAuthentication):
    def __init__(
        self,
        credential=None,
        resource_id="https://management.azure.com/.default",
        **kwargs,
    ):
        """Adapt any azure-identity credential to work with SDK that needs azure.common.credentials or msrestazure.
        Default resource is ARM (syntax of endpoint v2)
        :param credential: Any azure-identity credential (DefaultAzureCredential by default)
        :param str resource_id: The scope to use to get the token (default ARM)
        """
        try:
            from azure.core.pipeline.policies import BearerTokenCredentialPolicy
            from azure.identity import DefaultAzureCredential

        except ImportError:
            raise WorkflowError(
                "The Python 3 packages 'azure-core' and 'azure-identity' are required"
            )

        super(AzureIdentityCredentialAdapter, self).__init__(None)
        if credential is None:
            credential = DefaultAzureCredential()
        self._policy = BearerTokenCredentialPolicy(credential, resource_id, **kwargs)

    def _make_request(self):
        try:
            from azure.core.pipeline import PipelineContext, PipelineRequest
            from azure.core.pipeline.transport import HttpRequest
        except ImportError:
            raise WorkflowError("The Python 3 package azure-core is required")

        return PipelineRequest(
            HttpRequest("AzureIdentityCredentialAdapter", "https://fakeurl"),
            PipelineContext(None),
        )

    def set_token(self):
        """Ask the azure-core BearerTokenCredentialPolicy policy to get a token.
        Using the policy gives us for free the caching system of azure-core.
        We could make this code simpler by using private method, but by definition
        I can't assure they will be there forever, so mocking a fake call to the policy
        to extract the token, using 100% public API."""
        request = self._make_request()
        self._policy.on_request(request)
        # Read Authorization, and get the second part after Bearer
        token = request.http_request.headers["Authorization"].split(" ", 1)[1]
        self.token = {"access_token": token}

    def signed_session(self, session=None):
        self.set_token()
        return super(AzureIdentityCredentialAdapter, self).signed_session(session)


class AzBatchExecutor(RemoteExecutor):
    "Azure Batch Executor"

    def __init__(
        self,
        workflow: WorkflowExecutorInterface,
        dag: DAGExecutorInterface,
        stats: StatsExecutorInterface,
        logger: LoggerExecutorInterface,
        jobname="snakejob.{name}.{jobid}.sh",
        container_image=None,
        regions=None,
        location=None,
        cache=False,
        max_status_checks_per_second=1,
        az_batch_account_url=None,
        az_batch_enable_autoscale=False,
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
            from azure.batch import BatchServiceClient
            from azure.batch.batch_auth import SharedKeyCredentials
            from azure.identity import DefaultAzureCredential
            from azure.mgmt.batch import BatchManagementClient

            from snakemake.remote.AzBlob import AzureStorageHelper

        except ImportError:
            raise WorkflowError(
                "The Python 3 packages 'azure-batch', 'azure-mgmt-batch', and 'azure-identity'"
                " must be installed to use Azure Batch"
            )

        AZURE_BATCH_RESOURCE_ENDPOINT = "https://batch.core.windows.net/"

        # Here we validate that az blob credential is SAS
        # token because it is specific to azure batch executor
        self.validate_az_blob_credential_is_sas()
        self.azblob_helper = AzureStorageHelper()

        # get container from remote prefix
        self.prefix_container = str.split(workflow.default_remote_prefix, "/")[0]

        # setup batch configuration sets self.az_batch_config
        self.batch_config = AzBatchConfig(az_batch_account_url)
        logger.debug(f"AzBatchConfig: {self.mask_batch_config_as_string()}")

        self.workflow = workflow

        # handle case on OSX with /var/ symlinked to /private/var/ causing
        # issues with workdir not matching other workflow file dirs
        dirname = os.path.dirname(self.workflow.persistence.path)
        osxprefix = "/private"
        if osxprefix in dirname:
            dirname = dirname.removeprefix(osxprefix)

        self.workdir = dirname

        # Prepare workflow sources for build package
        self._set_workflow_sources()

        # Pool ids can only contain any combination of alphanumeric characters along with dash and underscore
        ts = datetime.datetime.now().strftime("%Y-%m%dT%H-%M-%S")
        self.pool_id = f"snakepool-{ts:s}"
        self.job_id = f"snakejob-{ts:s}"

        self.envvars = list(self.workflow.envvars) or []

        self.container_image = container_image or get_container_image()

        # enable autoscale flag
        self.az_batch_enable_autoscale = az_batch_enable_autoscale

        # Package workflow sources files and upload to storage
        self._build_packages = set()
        targz = self._generate_build_source_package()

        # removed after job failure/success
        self.resource_file = self._upload_build_source_package(
            targz, resource_prefix=self.batch_config.resource_file_prefix
        )

        # authenticate batch client from SharedKeyCredentials
        if (
            self.batch_config.batch_account_key is not None
            and self.batch_config.managed_identity_client_id is None
        ):
            logger.debug("Using batch account key for authentication...")
            creds = SharedKeyCredentials(
                self.batch_config.batch_account_name,
                self.batch_config.batch_account_key,
            )
        # else authenticate with managed indentity client id
        elif self.batch_config.managed_identity_client_id is not None:
            logger.debug("Using managed identity batch authentication...")
            creds = DefaultAzureCredential(
                managed_identity_client_id=self.batch_config.managed_identity_client_id
            )
            creds = AzureIdentityCredentialAdapter(
                credential=creds, resource_id=AZURE_BATCH_RESOURCE_ENDPOINT
            )

        self.batch_client = BatchServiceClient(
            creds, batch_url=self.batch_config.batch_account_url
        )

        if self.batch_config.managed_identity_resource_id is not None:
            self.batch_mgmt_client = BatchManagementClient(
                credential=DefaultAzureCredential(
                    managed_identity_client_id=self.batch_config.managed_identity_client_id
                ),
                subscription_id=self.batch_config.subscription_id,
            )

        try:
            self.create_batch_pool()
        except WorkflowError:
            logger.debug("Error: Failed to create batch pool, shutting down.")
            self.shutdown()

        try:
            self.create_batch_job()
        except WorkflowError:
            logger.debug("Error: Failed to create batch job, shutting down.")
            self.shutdown()

    def shutdown(self):
        # perform additional steps on shutdown
        # if necessary (jobs were cancelled already)

        logger.debug("Deleting AzBatch job")
        self.batch_client.job.delete(self.job_id)

        logger.debug("Deleting AzBatch pool")
        self.batch_client.pool.delete(self.pool_id)

        logger.debug("Deleting workflow sources from blob")

        self.azblob_helper.delete_from_container(
            self.prefix_container, self.resource_file.file_path
        )

        super().shutdown()

    def cancel(self):
        for task in self.batch_client.task.list(self.job_id):
            # strictly not need as job deletion also deletes task
            self.batch_client.task.terminate(self.job_id, task.id)
        self.shutdown()

    # mask_dict_vals masks sensitive keys from a dictionary of values for
    # logging used to mask dicts with sensitive information from logging
    @staticmethod
    def mask_dict_vals(mdict: dict, keys: list):
        ret_dict = mdict.copy()
        for k in keys:
            if k in ret_dict.keys() and ret_dict[k] is not None:
                ret_dict[k] = 10 * "*"
        return ret_dict

    # mask blob url is used to mask url values that may contain SAS
    # token information from being printed to the logs
    def mask_sas_urls(self, attrs: dict):
        attrs_new = attrs.copy()
        sas_pattern = r"\?[^=]+=([^?'\"]+)"
        mask = 10 * "*"

        for k, value in attrs.items():
            if value is not None and re.search(sas_pattern, str(value)):
                attrs_new[k] = re.sub(sas_pattern, mask, value)

        return attrs_new

    def mask_batch_config_as_string(self) -> str:
        masked_keys = self.mask_dict_vals(
            self.batch_config.__dict__,
            [
                "batch_account_key",
                "managed_identity_client_id",
            ],
        )
        masked_urls = self.mask_sas_urls(masked_keys)
        return pformat(masked_urls, indent=2)

    def run(
        self,
        job: ExecutorJobInterface,
        callback=None,
        submit_callback=None,
        error_callback=None,
    ):
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
        exec_job = f"/bin/bash -c 'tar xzf {self.resource_file.file_path} && {shlex.quote(exec_job)}'"

        # A string that uniquely identifies the Task within the Job.
        task_uuid = str(uuid.uuid1())
        task_id = f"{job.rule.name}-{task_uuid}"

        # This is the admin user who runs the command inside the container.
        user = batchmodels.AutoUserSpecification(
            scope=batchmodels.AutoUserScope.pool,
            elevation_level=batchmodels.ElevationLevel.admin,
        )

        # This is the docker image we want to run
        task_container_settings = batchmodels.TaskContainerSettings(
            image_name=self.container_image, container_run_options="--rm"
        )

        # https://docs.microsoft.com/en-us/python/api/azure-batch/azure.batch.models.taskaddparameter?view=azure-python
        # all directories recursively below the AZ_BATCH_NODE_ROOT_DIR (the root of Azure Batch directories on the node)
        # are mapped into the container, all Task environment variables are mapped into the container,
        # and the Task command line is executed in the container
        task = batch.models.TaskAddParameter(
            id=task_id,
            command_line=exec_job,
            container_settings=task_container_settings,
            resource_files=[self.resource_file],  # Snakefile, envs, yml files etc.
            user_identity=batchmodels.UserIdentity(auto_user=user),
            environment_settings=envsettings,
        )

        # register job as active, using your own namedtuple.
        self.batch_client.task.add(self.job_id, task)
        self.active_jobs.append(
            AzBatchJob(job, self.job_id, task_id, callback, error_callback)
        )
        logger.debug(f"Added AzBatch task {task_id}")
        logger.debug(
            f"Added AzBatch task {pformat(self.mask_sas_urls(task.__dict__), indent=2)}"
        )

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
        except Exception:
            content = ""

        return content

    async def _wait_for_jobs(self):
        import azure.batch.models as batchmodels

        while True:
            # always use self.lock to avoid race conditions
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()

            # Loop through active jobs and act on status
            for batch_job in active_jobs:
                async with self.status_rate_limiter:
                    logger.debug(f"Monitoring {len(active_jobs)} active AzBatch tasks")
                    task = self.batch_client.task.get(self.job_id, batch_job.task_id)

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
                        logger.debug(
                            "task {} completed: result={} exit_code={}\n".format(
                                batch_job.task_id, task.execution_info.result, rc
                            )
                        )
                        logger.debug(
                            "task {} completed: run_time={}, retry_count={}\n".format(
                                batch_job.task_id, str(dt), rt
                            )
                        )

                        def print_output():
                            logger.debug(
                                "task {}: stderr='{}'\n".format(
                                    batch_job.task_id, stderr
                                )
                            )
                            logger.debug(
                                "task {}: stdout='{}'\n".format(
                                    batch_job.task_id, stdout
                                )
                            )

                        if (
                            task.execution_info.result
                            == batchmodels.TaskExecutionResult.failure
                        ):
                            logger.error(
                                f"Azure task failed: code={str(task.execution_info.failure_info.code)}, message={str(task.execution_info.failure_info.message)}"
                            )
                            for d in task.execution_info.failure_info.details:
                                logger.error(f"Error Details: {str(d)}")
                            print_output()
                            batch_job.error_callback(batch_job.job)
                        elif (
                            task.execution_info.result
                            == batchmodels.TaskExecutionResult.success
                        ):
                            batch_job.callback(batch_job.job)
                        else:
                            logger.error(
                                "Unknown Azure task execution result: {}".format(
                                    task.execution_info.result
                                )
                            )
                            print_output()
                            batch_job.error_callback(batch_job.job)

                    # The operation is still running
                    else:
                        logger.debug(
                            f"task {batch_job.task_id}: creation_time={task.creation_time} state={task.state} node_info={task.node_info}\n"
                        )
                        still_running.append(batch_job)

                        # fail if start task fails on a node or node state becomes unusable
                        # and stream stderr stdout to stream
                        node_list = self.batch_client.compute_node.list(self.pool_id)
                        for n in node_list:
                            # error on unusable node (this occurs if your container image fails to pull)
                            if n.state == "unusable":
                                if n.errors is not None:
                                    for e in n.errors:
                                        logger.error(
                                            f"Azure task error: {e.message}, {e.error_details[0].__dict__}"
                                        )
                                logger.error(
                                    "A node entered an unusable state, quitting."
                                )
                                return

                            if n.start_task_info is not None and (
                                n.start_task_info.result
                                == batchmodels.TaskExecutionResult.failure
                            ):
                                try:
                                    stderr_file = (
                                        self.batch_client.file.get_from_compute_node(
                                            self.pool_id, n.id, "/startup/stderr.txt"
                                        )
                                    )
                                    stderr_stream = self._read_stream_as_string(
                                        stderr_file, "utf-8"
                                    )
                                except Exception:
                                    stderr_stream = ""

                                try:
                                    stdout_file = (
                                        self.batch_client.file.get_from_compute_node(
                                            self.pool_id, n.id, "/startup/stdout.txt"
                                        )
                                    )
                                    stdout_stream = self._read_stream_as_string(
                                        stdout_file, "utf-8"
                                    )
                                except Exception:
                                    stdout_stream = ""

                                logger.error(
                                    "Azure start task execution failed on node: {}.\nSTART_TASK_STDERR:{}\nSTART_TASK_STDOUT: {}".format(
                                        n.start_task_info.failure_info.message,
                                        stdout_stream,
                                        stderr_stream,
                                    )
                                )
                                return

            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            await sleep()

    def create_batch_pool(self):
        """Creates a pool of compute nodes"""

        import azure.batch._batch_service_client as bsc
        import azure.batch.models as batchmodels
        import azure.mgmt.batch.models as mgmtbatchmodels

        image_ref = bsc.models.ImageReference(
            publisher=self.batch_config.batch_pool_image_publisher,
            offer=self.batch_config.batch_pool_image_offer,
            sku=self.batch_config.batch_pool_image_sku,
            version="latest",
        )

        # optional subnet network configuration
        # requires AAD batch auth insead of batch key auth
        network_config = None
        if self.batch_config.batch_pool_subnet_id is not None:
            network_config = batchmodels.NetworkConfiguration(
                subnet_id=self.batch_config.batch_pool_subnet_id
            )

        # configure a container registry

        # Specify container configuration, fetching an image
        #  https://docs.microsoft.com/en-us/azure/batch/batch-docker-container-workloads#prefetch-images-for-container-configuration
        container_config = batchmodels.ContainerConfiguration(
            type="dockerCompatible", container_image_names=[self.container_image]
        )

        user = None
        passw = None
        identity_ref = None
        registry_conf = None

        if self.batch_config.container_registry_url is not None:
            if (
                self.batch_config.container_registry_user is not None
                and self.batch_config.container_registry_pass is not None
            ):
                user = self.batch_config.container_registry_user
                passw = self.batch_config.container_registry_pass
            elif self.batch_config.managed_identity_resource_id is not None:
                identity_ref = batchmodels.ComputeNodeIdentityReference(
                    resource_id=self.batch_config.managed_identity_resource_id
                )
            else:
                raise WorkflowError(
                    "No container registry authentication scheme set. Please set the BATCH_CONTAINER_REGISTRY_USER and BATCH_CONTAINER_REGISTRY_PASS or set MANAGED_IDENTITY_CLIENT_ID and MANAGED_IDENTITY_RESOURCE_ID."
                )

            registry_conf = [
                batchmodels.ContainerRegistry(
                    registry_server=self.batch_config.container_registry_url,
                    identity_reference=identity_ref,
                    user_name=str(user),
                    password=str(passw),
                )
            ]

            # Specify container configuration, fetching an image
            #  https://docs.microsoft.com/en-us/azure/batch/batch-docker-container-workloads#prefetch-images-for-container-configuration
            container_config = batchmodels.ContainerConfiguration(
                type="dockerCompatible",
                container_image_names=[self.container_image],
                container_registries=registry_conf,
            )

        # default to no start task
        start_task = None

        # if configured us start task bash script from sas url
        if self.batch_config.batch_node_start_task_sasurl is not None:
            _SIMPLE_TASK_NAME = "start_task.sh"
            start_task_admin = batchmodels.UserIdentity(
                auto_user=batchmodels.AutoUserSpecification(
                    elevation_level=batchmodels.ElevationLevel.admin,
                    scope=batchmodels.AutoUserScope.pool,
                )
            )
            start_task = batchmodels.StartTask(
                command_line=f"bash {_SIMPLE_TASK_NAME}",
                resource_files=[
                    batchmodels.ResourceFile(
                        file_path=_SIMPLE_TASK_NAME,
                        http_url=self.batch_config.batch_node_start_task_sasurl,
                    )
                ],
                user_identity=start_task_admin,
            )

        # autoscale requires the initial dedicated node count to be zero
        if self.az_batch_enable_autoscale:
            self.batch_config.batch_pool_node_count = 0

        node_communication_strategy = None
        if self.batch_config.batch_node_communication_mode is not None:
            node_communication_strategy = batchmodels.NodeCommunicationMode.simplified

        new_pool = batchmodels.PoolAddParameter(
            id=self.pool_id,
            virtual_machine_configuration=batchmodels.VirtualMachineConfiguration(
                image_reference=image_ref,
                container_configuration=container_config,
                node_agent_sku_id=self.batch_config.batch_pool_vm_node_agent_sku_id,
            ),
            network_configuration=network_config,
            vm_size=self.batch_config.batch_pool_vm_size,
            target_dedicated_nodes=self.batch_config.batch_pool_node_count,
            target_node_communication_mode=node_communication_strategy,
            target_low_priority_nodes=0,
            start_task=start_task,
            task_slots_per_node=self.batch_config.batch_tasks_per_node,
            task_scheduling_policy=batchmodels.TaskSchedulingPolicy(
                node_fill_type=self.batch_config.batch_node_fill_type
            ),
        )

        # create pool if not exists
        try:
            logger.debug(f"Creating pool: {self.pool_id}")
            self.batch_client.pool.add(new_pool)

            if self.az_batch_enable_autoscale:
                # define the autoscale formula
                formula = """$samples = $PendingTasks.GetSamplePercent(TimeInterval_Minute * 5);
                            $tasks = $samples < 70 ? max(0,$PendingTasks.GetSample(1)) : max( $PendingTasks.GetSample(1), avg($PendingTasks.GetSample(TimeInterval_Minute * 5)));
                            $targetVMs = $tasks > 0? $tasks:max(0, $TargetDedicatedNodes/2);
                            $TargetDedicatedNodes = max(0, min($targetVMs, 10));
                            $NodeDeallocationOption = taskcompletion;"""

                # Enable autoscale; specify the formula
                self.batch_client.pool.enable_auto_scale(
                    self.pool_id,
                    auto_scale_formula=formula,
                    # the minimum allowed autoscale interval is 5 minutes
                    auto_scale_evaluation_interval=datetime.timedelta(minutes=5),
                    pool_enable_auto_scale_options=None,
                    custom_headers=None,
                    raw=False,
                )

            # update pool with managed identity, enables batch nodes to act as managed identity
            if self.batch_config.managed_identity_resource_id is not None:
                mid = mgmtbatchmodels.BatchPoolIdentity(
                    type=mgmtbatchmodels.PoolIdentityType.user_assigned,
                    user_assigned_identities={
                        self.batch_config.managed_identity_resource_id: mgmtbatchmodels.UserAssignedIdentities()
                    },
                )
                params = mgmtbatchmodels.Pool(identity=mid)
                self.batch_mgmt_client.pool.update(
                    resource_group_name=self.batch_config.resource_group,
                    account_name=self.batch_config.batch_account_name,
                    pool_name=self.pool_id,
                    parameters=params,
                )

        except batchmodels.BatchErrorException as err:
            if err.error.code != "PoolExists":
                raise WorkflowError(
                    f"Error: Failed to create pool: {err.error.message}"
                )
            else:
                logger.debug(f"Pool {self.pool_id} exists.")

    def create_batch_job(self):
        """Creates a job with the specified ID, associated with the specified pool"""
        import azure.batch._batch_service_client as bsc

        logger.debug(f"Creating job {self.job_id}")

        self.batch_client.job.add(
            bsc.models.JobAddParameter(
                id=self.job_id,
                constraints=bsc.models.JobConstraints(max_task_retry_count=0),
                pool_info=bsc.models.PoolInformation(pool_id=self.pool_id),
            )
        )

    @staticmethod
    def validate_az_blob_credential_is_sas():
        """
        Validates that the AZ_BLOB_CREDENTIAL is a valid storage account SAS
        token, required when using --az-batch with AzBlob remote.
        """
        cred = os.environ.get("AZ_BLOB_CREDENTIAL")
        if cred is not None:
            # regex pattern for Storage Account SAS Token
            rgx = r"\?sv=.*&ss=.*&srt=.*&sp=.*&se=.*&st=.*&spr=.*&sig=.*"
            if re.compile(rgx).match(cred) is None:
                raise WorkflowError(
                    "AZ_BLOB_CREDENTIAL is not a valid storage account SAS token."
                )

    # from google_lifesciences.py
    def _set_workflow_sources(self):
        """We only add files from the working directory that are config related
        (e.g., the Snakefile or a config.yml equivalent), or checked into git.
        """
        self.workflow_sources = []

        for wfs in self.dag.get_sources():
            if os.path.isdir(wfs):
                for dirpath, dirnames, filenames in os.walk(wfs):
                    self.workflow_sources.extend(
                        [check_source_size(os.path.join(dirpath, f)) for f in filenames]
                    )
            else:
                self.workflow_sources.append(check_source_size(os.path.abspath(wfs)))

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
        targz = os.path.join(tempfile.gettempdir(), f"snakemake-{tmpname}.tar.gz")
        tar = tarfile.open(targz, "w:gz")

        # Add all workflow_sources files
        for filename in self.workflow_sources:
            arcname = filename.replace(self.workdir + os.path.sep, "")
            tar.add(filename, arcname=arcname)
        logger.debug(
            f"Created {targz} with the following contents: {self.workflow_sources}"
        )
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

    def _upload_build_source_package(self, targz, resource_prefix=""):
        """given a .tar.gz created for a workflow, upload it to the blob
        storage account, only if the blob doesn't already exist.
        """

        import azure.batch.models as batchmodels

        blob_name = os.path.join(resource_prefix, os.path.basename(targz))

        # upload blob to storage using storage helper
        bc = self.azblob_helper.upload_to_azure_storage(
            self.prefix_container, targz, blob_name=blob_name
        )

        # return resource file
        return batchmodels.ResourceFile(http_url=bc.url, file_path=blob_name)
