__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import shlex
import sys
from collections import namedtuple

from snakemake.executors import ClusterExecutor, sleep
from snakemake.executors.common import format_cli_arg, join_cli_args
from snakemake.logging import logger
from snakemake.resources import DefaultResources
from snakemake.common import async_lock

# Just import flux once
try:
    import flux
    import flux.job
    from flux.job import JobspecV1
except ImportError:
    flux = None

FluxJob = namedtuple(
    "FluxJob", "job jobname jobid callback error_callback flux_future flux_logfile"
)


class FluxExecutor(ClusterExecutor):
    """
    The Flux executor deploys workflows to a flux cluster.
    """

    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob.{name}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
    ):
        super().__init__(
            workflow,
            dag,
            None,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            assume_shared_fs=False,
            max_status_checks_per_second=10,
        )

        # Attach variables for easy access
        self.workdir = os.path.realpath(os.path.dirname(self.workflow.persistence.path))
        self.envvars = list(self.workflow.envvars) or []

        # Quit early if we can't access the flux api
        if not flux:
            raise WorkflowError(
                "Cannot import flux. Is it installed (https://flux-framework.org) and available to you with Python bindings?"
            )
        self._fexecutor = flux.job.FluxExecutor()

    def cancel(self):
        """
        cancel execution, usually by way of control+c. Cleanup is done in
        shutdown (deleting cached workdirs in Google Cloud Storage
        """
        for job in self.active_jobs:
            if not job.flux_future.done():
                flux.job.cancel(self.f, job.jobid)
        self.shutdown()

    def _set_job_resources(self, job):
        """
        Given a particular job, generate the resources that it needs,
        including default regions and the virtual machine configuration
        """
        self.default_resources = DefaultResources(
            from_other=self.workflow.default_resources
        )

    def get_snakefile(self):
        assert os.path.exists(self.workflow.main_snakefile)
        return self.workflow.main_snakefile

    def _get_jobname(self, job):
        # Use a dummy job name (human readable and also namespaced)
        return "snakejob-%s-%s-%s" % (self.run_namespace, job.name, job.jobid)

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        """
        Submit a job to flux.
        """
        super()._run(job)

        flux_logfile = job.logfile_suggestion(".snakemake/flux_logs")
        os.makedirs(os.path.dirname(flux_logfile), exist_ok=True)

        # Prepare job resourcces
        self._set_job_resources(job)

        # The entire snakemake command to run, etc
        command = self.format_job_exec(job)
        logger.debug(command)

        # Generate the flux job
        # flux does not support mem_mb, disk_mb
        fluxjob = JobspecV1.from_command(command=shlex.split(command))

        # A duration of zero (the default) means unlimited
        fluxjob.duration = job.resources.get("runtime", 0)
        fluxjob.stderr = flux_logfile

        # Ensure the cwd is the snakemake working directory
        fluxjob.cwd = self.workdir
        fluxjob.environment = dict(os.environ)
        flux_future = self._fexecutor.submit(fluxjob)

        logger.info(
            f"Job {job.jobid} has been submitted with flux jobid {flux_future.jobid()} (log: {flux_logfile})."
        )

        # Waiting for the jobid is a small performance penalty, same as calling flux.job.submit
        self.active_jobs.append(
            FluxJob(
                job,
                str(job),
                str(flux_future.jobid()),
                callback,
                error_callback,
                flux_future,
                flux_logfile,
            )
        )

    async def _wait_for_jobs(self):
        """
        Wait for jobs to complete. This means requesting their status,
        and then marking them as finished when a "done" parameter
        shows up. Even for finished jobs, the status should still return
        """
        while True:
            # always use self.lock to avoid race conditions
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()

            # Loop through active jobs and act on status
            for j in active_jobs:
                logger.debug("Checking status for job {}".format(j.jobid))
                if j.flux_future.done():
                    # The exit code can help us determine if the job was successful
                    try:
                        exit_code = j.flux_future.result(0)
                    except RuntimeError:
                        # job did not complete
                        self.print_job_error(j.job, jobid=j.jobid)
                        j.error_callback(j.job)

                    else:
                        # the job finished (but possibly with nonzero exit code)
                        if exit_code != 0:
                            self.print_job_error(
                                j.job,
                                jobid=j.jobid,
                                aux_logs=[j.flux_logfile],
                            )
                            j.error_callback(j.job)
                            continue

                        # Finished and success!
                        j.callback(j.job)

                # Otherwise, we are still running
                else:
                    still_running.append(j)
            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)

            # Sleeps for 10 seconds
            await sleep()
