__author__ = "Andreas Wilm Köster"
__copyright__ = "Copyright 2020, Andreas Wilm"
__email__ = "andreas.wilm@gmail.com"
__license__ = "MIT"

import os
from collections import namedtuple
import datetime
import uuid

from snakemake.executors import ClusterExecutor, sleep
from snakemake.exceptions import WorkflowError

# in azure batch parlance this is actually a task
# FIXME should probably make job_id part of this to get rid of self.job_id
AzBatchJob = namedtuple(
    "AzBatchJob", "job task_id callback error_callback"
)


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
        latency_wait=3,
        local_input=None,
        restart_times=None,
        exec_job=None,
        max_status_checks_per_second=1,
        az_batch_config=None):


        # FIXME from google_lifesciences.py
        # Attach variables for easy access
        self.workflow = workflow
        #self.quiet = quiet
        self.workdir = os.path.dirname(self.workflow.persistence.path)
        #self._save_storage_cache = cache

        # Pool ids can only contain any combination of alphanumeric characters along with dash and underscore
        ts = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        self.pool_id = "snakemakepool-{:s}".format(ts)
        self.job_id = "snakemakejob-{:s}".format(ts)

        # Relative path for running on instance
        self._set_snakefile()

        # Prepare workflow sources for build package
        # FIXME google_lifesciences.py: self._set_workflow_sources()

        exec_job = exec_job or (
            "snakemake {target} --snakefile %s "
            "--force -j{cores} --keep-target-files --keep-remote "
            "--latency-wait 0 "
            "--attempt 1 {use_threads} "
            "{overwrite_config} {rules} --nocolor "
            "--notemp --no-hooks --nolock " % self.snakefile
        )

        try:
            import azure.batch._batch_service_client as batch# https://docs.microsoft.com/en-us/azure/batch/quick-run-python uses old 6.0 API without _
            import azure.batch.batch_auth as batch_auth
            #import azure.batch.models as batchmodels
        except ImportError:
            raise WorkflowError(
                "The Python 3 package 'azure-batch'"
                " must be installed to use Azure Batch"
            )

        self.az_batch_config = az_batch_config

        # Create a Batch service client.
        credentials = batch_auth.SharedKeyCredentials(
            self.az_batch_config['BATCH_ACCOUNT_NAME'],
            self.az_batch_config['BATCH_ACCOUNT_KEY'])
        self.batch_client = batch.BatchServiceClient(
            credentials, batch_url=self.az_batch_config['BATCH_ACCOUNT_URL'])

        # Create the pool that will contain the compute nodes that will execute the
        # tasks.
        print("FIXME creating pool")    
        self.create_pool(self.batch_client, self.pool_id)

        # Create the job that will run the tasks.
        print("FIXME creating job")    
        self.create_job(self.batch_client, self.job_id, self.pool_id)
        
        super().__init__(
            workflow,
            dag,
            None,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            restart_times=restart_times,
            exec_job=exec_job,
            assume_shared_fs=False,
            max_status_checks_per_second=10,
        )


    def shutdown(self):
        # perform additional steps on shutdown if necessary (jobs were cancelled already) 
        print("FIXME deleting job")
        self.batch_client.job.delete(self.job_id)
        print("FIXME deleting pool")
        self.batch_client.pool.delete(self.pool_id)
        super().shutdown()
        # FIXME from google_lifesciences.py


    def cancel(self):
        for task in self.batch_client.task.list(self.job_id):
            print("FIXME cancel active jobs here")
        self.shutdown()


    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):

        import azure.batch._batch_service_client as batch# https://docs.microsoft.com/en-us/azure/batch/quick-run-python uses old 6.0 API without _

        super()._run(job)

        # obtain job execution command
        use_threads = "--force-use-threads" if not job.is_group() else ""
        exec_job = self.format_job(
            self.exec_job, job, _quote_all=True,
            use_threads=use_threads)

        # the job is called task in Azure Batch parlance
        # whereas a job is a group of tasks     
        task_id = str(uuid.uuid1())# A string that uniquely identifies the Task within the Job. 
        task = batch.models.TaskAddParameter(id=task_id, command_line=exec_job)
        # FIXME need resource_files?
        # FIXME need environment_settings as well?

        # register job as active, using your own namedtuple.
        # The namedtuple must at least contain the attributes
        # job, jobid, callback, error_callback.
        #taskname = "FIXME"
        #self.active_jobs.append(
        #    AzBatchJob(task, task_id, taskname, callback, error_callback))
        # FIXME could we get rid of AzBatchJob altogether in favour of self.batch_service_client?
        # looks like we need to use batch_service_client to actually run anything 
        # FIXME autorestart should be handled here
        self.batch_client.task.add(self.job_id, task)
        self.active_jobs.append(
            AzBatchJob(job, task_id, callback, error_callback)
        )
        print("FIXME added task with id {}".format(task_id))


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

            print("FIXME {:d} tasks in list".format(len(active_jobs)))
            # Loop through active jobs and act on status
            for batch_job in active_jobs:
                print("FIXME checking task {}".format(batch_job.task_id))
                # use self.status_rate_limiter to avoid too many API calls.
                with self.status_rate_limiter:
                    task = self.batch_client.task.get(self.job_id, batch_job.task_id)
                    #if j.task.status == batchmodels.TaskState.completed:
                    if task.state == batchmodels.TaskState.completed:
                        import pdb; pdb.set_trace()
                        print("FIXME ended task {} has result {}".format(batch_job.task_id, task.execution_info.result))
                        if task.execution_info.result == batchmodels.TaskExecutionResult.failure:
                            batch_job.error_callback(batch_job.job)
                        elif task.execution_info.result == batchmodels.TaskExecutionResult.success:
                            batch_job.callback(batch_job.job)
                        else:
                            raise ValueError("Unknown task execution result: {}".format(task.execution_info.result))

                    # The operation is still running
                    else:
                        print("FIXME {} still running".format(batch_job.task_id))
                        still_running.append(batch_job)

            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


    def create_pool(self, batch_service_client, pool_id):
        """Creates a pool of compute nodes 
        """
        
        import azure.batch.models as batchmodels
        import azure.batch._batch_service_client as batch# https://docs.microsoft.com/en-us/azure/batch/quick-run-python uses old 6.0 API without _

        print('DEBUG Creating pool [{}]...'.format(pool_id))# FIXME AW

        # Create a new pool of Linux compute nodes using an Azure Virtual Machines
        # Marketplace image. For more information about creating pools of Linux
        # nodes, see:
        # https://azure.microsoft.com/documentation/articles/batch-linux-nodes/
        print("FIXME hardcoded image reference")
        print("FIXME  make kubernetes/singularity mandatory?")

        new_pool = batch.models.PoolAddParameter(
            id=pool_id,
            virtual_machine_configuration=batchmodels.VirtualMachineConfiguration(
                image_reference=batchmodels.ImageReference(
                    publisher="Canonical",
                    offer="UbuntuServer",
                    sku="18.04-LTS",
                    version="latest"
                    # FIXME make config
                ),
                node_agent_sku_id="batch.node.ubuntu 18.04"),
            # FIXME autoscaling
            vm_size = self.az_batch_config['BATCH_POOL_VM_SIZE'],
            target_dedicated_nodes = self.az_batch_config['BATCH_POOL_NODE_COUNT']
        )
        batch_service_client.pool.add(new_pool)


    @staticmethod
    def create_job(batch_service_client, job_id, pool_id):
        """Creates a job with the specified ID, associated with the specified pool
        """
        import azure.batch._batch_service_client as batch# https://docs.microsoft.com/en-us/azure/batch/quick-run-python uses old 6.0 API without _
        
        print('DEBUG Creating job [{}]...'.format(job_id))# FIXME AW

        job = batch.models.JobAddParameter(
            id=job_id,
            pool_info=batch.models.PoolInformation(pool_id=pool_id))

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
