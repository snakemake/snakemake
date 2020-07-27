__author__ = "Andreas Wilm Köster"
__copyright__ = "Copyright 2020, Andreas Wilm"
__email__ = "andreas.wilm@gmail.com"
__license__ = "MIT"

import os
from collections import namedtuple

from snakemake.executors import ClusterExecutor, sleep
from snakemake.exceptions import WorkflowError


AzBatchJob = namedtuple(
    "AzBatchJob", "job jobname jobid callback error_callback"
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
            credentials,
            batch_url=self.az_batch_config['BATCH_ACCOUNT_URL'])

        # Create the pool that will contain the compute nodes that will execute the
        # tasks.
        print("FIXME creating pool")    
        self.create_pool(self.batch_client, self.az_batch_config['BATCH_POOL_ID'])

        # Create the job that will run the tasks.
        print("FIXME creating job")    
        self.create_job(self.batch_client, self.az_batch_config['BATCH_JOB_ID'], 
            self.az_batch_config['BATCH_POOL_ID'])
        
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
        self.batch_client.job.delete(self.az_batch_config['BATCH_JOB_ID'])
        print("FIXME deleting pool")
        self.batch_client.pool.delete(self.az_batch_config['BATCH_POOL_ID'])
        super().shutdown()


    def cancel(self):
        for job in self.active_jobs:
            print("FIXME cancel active jobs here")
        self.shutdown()


    def run(self, job,
            callback=None,
            submit_callback=None,
            error_callback=None):

        super()._run(job)
        # obtain job execution command
        exec_job = self.format_job(
            self.exec_job, job, _quote_all=True,
            use_threads="--force-use-threads" if not job.is_group() else "")

        print("FIXME submit job here, and obtain job ids from the backend")

        # register job as active, using your own namedtuple.
        # The namedtuple must at least contain the attributes
        # job, jobid, callback, error_callback.
        # FIXME self.active_jobs.append(AzBatchJob(
        #    job, jobid, callback, error_callback))


    def _wait_for_jobs(self):
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
                with self.status_rate_limiter:
                    pass
                    # Retrieve status of job j from your backend via j.jobid
                    # Handle completion and errors, calling either j.callback(j.job)
                    # or j.error_callback(j.job)
                    # In case of error, add job j to still_running.
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()


    def create_pool(self, batch_service_client, pool_id):
        """
        Creates a pool of compute nodes with the specified OS settings.
        :param batch_service_client: A Batch service client.
        :type batch_service_client: `azure.batch.BatchServiceClient`
        :param str pool_id: An ID for the new pool.
        :param str publisher: Marketplace image publisher
        :param str offer: Marketplace image offer
        :param str sku: Marketplace image sku
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
                ),
                node_agent_sku_id="batch.node.ubuntu 18.04"),
            # FIXME autoscaling
            vm_size=self.az_batch_config['BATCH_POOL_VM_SIZE'],
            target_dedicated_nodes=self.az_batch_config['BATCH_POOL_NODE_COUNT']
        )
        batch_service_client.pool.add(new_pool)


    @staticmethod
    def create_job(batch_service_client, job_id, pool_id):
        """
        Creates a job with the specified ID, associated with the specified pool.
        :param batch_service_client: A Batch service client.
        :type batch_service_client: `azure.batch.BatchServiceClient`
        :param str job_id: The ID for the job.
        :param str pool_id: The ID for the pool.
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
