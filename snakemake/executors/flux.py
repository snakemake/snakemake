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

# Just import flux once
try:
    import flux
    import flux.job
except ImportError:
    flux = None

FluxJob = namedtuple("FluxJob", "job jobname jobid callback error_callback flux_future")


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
        self._init_flux()

    def _prepare_job_formatter(self):
        """
        A job doesn't have an easy status command - instead we do a job
        listing with a particular format. See src/cmd/flux-jobs.py#L44
        in flux-framework/flux-core for more attributes.
        """
        jobs_format = (
            "{id.f58:>12} {username:<8.8} {name:<10.10} {status_abbrev:>2.2} "
            "{ntasks:>6} {nnodes:>6h} {runtime!F:>8h} {success} {exception.occurred}"
            "{exception.note} {exception.type} {result} {runtime} {status}"
            "{ranks:h} {t_remaining} {annotations}"
        )
        self.jobs_formatter = flux.job.JobInfoFormat(jobs_format)

        # Note there is no attr for "id", its always returned
        fields2attrs = {
            "id.f58": (),
            "username": ("userid",),
            "exception.occurred": ("exception_occurred",),
            "exception.type": ("exception_type",),
            "exception.note": ("exception_note",),
            "runtime": ("t_run", "t_cleanup"),
            "status": ("state", "result"),
            "status_abbrev": ("state", "result"),
            "t_remaining": ("expiration", "state", "result"),
        }

        # Set job attributes we will use later to get job statuses
        self.job_attrs = set()
        for field in self.jobs_formatter.fields:
            if field not in fields2attrs:
                self.job_attrs.update((field,))
            else:
                self.job_attrs.update(fields2attrs[field])

    def _init_flux(self):
        """
        Ensure we can import flux in the environment.
        """
        if not flux:
            sys.exit(
                "Cannot import flux. Is a cluster available to you with Python bindings?"
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
        return self.workflow.main_snakefile.replace(self.workdir, "").strip(os.sep)

    def _get_jobname(self, job):
        # Use a dummy job name (human readable and also namespaced)
        return "snakejob-%s-%s-%s" % (self.run_namespace, job.name, job.jobid)

    def _generate_environment(self):
        """
        Loop through envvars (keys to host environment) and add
        any that are requested for the runtime environment
        """
        envvars = dict(os.environ)
        for key in self.envvars:
            try:
                envvars[key] = os.environ[key]
            except Exception:
                continue
        return envvars

    def format_job_exec(self, job):
        """
        We don't want a python -m to run.
        """
        prefix = self.get_job_exec_prefix(job)
        if prefix:
            prefix += " &&"
        suffix = self.get_job_exec_suffix(job)
        if suffix:
            suffix = f"&& {suffix}"
        return join_cli_args(
            [
                prefix,
                self.get_envvar_declarations(),
                "snakemake",
                format_cli_arg("--snakefile", self.get_snakefile()),
                self.get_job_args(job),
                self.general_args,
                suffix,
            ]
        )

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        """
        Submit a job to flux.
        """
        super()._run(job)

        from flux.job import JobspecV1

        # This needs to be called after the line above
        self._prepare_job_formatter()

        # Prepare job resourcces
        self._set_job_resources(job)

        # The entire snakemake command to run, etc
        command = self.format_job_exec(job)
        logger.debug(command)

        # Generate the flux job
        fluxjob = JobspecV1.from_command(command=shlex.split(command))

        # Ensure the cwd is the snakemake working directory
        fluxjob.cwd = self.workdir
        fluxjob.environment = self._generate_environment()
        flux_future = self._fexecutor.submit(fluxjob)

        # Waiting for the jobid is a small performance penalty, same as calling flux.job.submit
        self.active_jobs.append(
            FluxJob(
                job,
                str(job),
                str(flux_future.jobid()),
                callback,
                error_callback,
                flux_future,
            )
        )

    def _wait_for_jobs(self):
        """
        Wait for jobs to complete. This means requesting their status,
        and then marking them as finished when a "done" parameter
        shows up. Even for finished jobs, the status should still return
        """
        while True:
            # always use self.lock to avoid race conditions
            with self.lock:
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
                            self.print_job_error(j.job, jobid=j.jobid)
                            j.error_callback(j.job)
                            continue

                        # Finished and success!
                        j.callback(j.job)

                # Otherwise, we are still running
                else:
                    still_running.append(j)
            with self.lock:
                self.active_jobs.extend(still_running)

            # Sleeps for 10 seconds
            sleep()
