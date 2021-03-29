import os
import re
import sys
import time
import shutil
import tarfile
import tempfile

from snakemake.jobs import Job
from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.exceptions import print_exception
from snakemake.exceptions import log_verbose_traceback
from snakemake.exceptions import WorkflowError
from snakemake.executors import ClusterExecutor
from snakemake.utils import makedirs
from snakemake.io import get_wildcard_names, Wildcards

SlurmJob = namedtuple("SlurmJob", "job jobid callback error_callback")


class SlurmExecutor(ClusterExecutor):
    """
    the SLURM_Executor abstracts execution on SLURM
    clusters using snakemake resource string
    """

    default_jobscript = "slurm_jobscript.py"

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
        container_image=None,
        regions=None,
        location=None,
        cache=False,
        local_input=None,
        restart_times=None,
        exec_job=None,
        max_status_checks_per_second=60,
    ):

        super().__init__(
            workflow,
            dag,
            cores,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            cluster_config=cluster_config,
            restart_times=restart_times,
            assume_shared_fs=assume_shared_fs,
            max_status_checks_per_second=max_status_checks_per_second,
        )
        self.exec_job += " --slurm-jobstep"

    def shutdown(self):
        """
        waiting for all jobs to be ready
        """
        super().shutdown()
        # perform additional steps on shutdown if necessary

    def cancel(self):
        for job in self.active_jobs:
            try:
                # timeout set to 60, because a scheduler cycle usually is
                # about 30 sec, but can be longer in extreme cases.
                # Under 'normal' circumstances, 'scancel' is executed in
                # virtually no time.
                subprocess.run('scancel {job.jobid}', timeout = 60)
            except TimeoutExpired:
                pass # shell we ignore the timeout, here?
        self.shutdown()

    def run_jobs(self, jobs, callback=None, submit_callback=None, error_callback=None):
        """Run a list of jobs that is ready at a given point in time.

        By default, this method just runs each job individually.
        This behavior is inherited and therefore this method can be removed from the skeleton if the
        default behavior is intended.
        This method can be overwritten to submit many jobs in a more efficient way than one-by-one.

        Note that in any case, for each job, the callback functions have to be called individually!
        """
        # if groupjob:
        #    attributes of groubjobs (snakemake/jobs.py class GroupJob , _resources
        #    - see topology

        for job in jobs:
            # TODO replace by own code
            self.run(
                job,
                callback=callback,  # to be executed upon ready jobs
                submit_callback=submit_callback,  # to be hold for --immediate flag - ignored within this executor
                error_callback=error_callback,  #
            )

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        # generic part of a submission string:
        call = "sbatch -A {account} -p {partition} -t {walltime_minutes} \
                -J {job.name} -o .snakemake/slurm_logs/%x_%j.log \
                --export=ALL"
        if job.resources.get("constraint"):
            call += " -C {constraint}"  # from (default)resources
        # TODO: implement when tempfs-resource is defined
        # if job.resources.get()
        if job.resources.get("mem_mb_per_cpu"):
            call += " --mem-per-cpu={mem_mb_per_cpu}"
        else:
            call += " --mem {mem_mb}"

        # if ressources == MPI:
        if job.resources.get("mpi", False):
            # submit MPI job
            # mold submissing string, submit according to generic cluster exec
            pass
        # ordinary smp application:
        elif not job.resources.get("mpi") and not job.is_group():
            call = " -n 1 -c {threads} {exec_job}".format(self.exec_job)
            try:
                out = subprocess.check_output(call)
            except:
                pass  # check template

            jobid = out.split(" ")[-1]
            self.active_jobs.append(SlurmJob(job, jobid, callback, error_callback))
        else: # job.is_group:
            pass

    def job_status(self):
        """
        obtain SLURM job status of submitted jobs
        """
        STATUS_ATTEMPS = 10
        for i in range(STATUS_ATTEMPS):
            # use self.status_rate_limiter to avoid too many API calls.
            with self.status_rate_limiter:
                try:
                    sacct_res = subprocess.check_output(
                        "sacct -P -b -j {} -n".format(jobid), shell=True
                    )
                    res = {
                        x.split("|")[0]: x.split("|")[1]
                        for x in sacct_res.decode().strip().split("\n")
                    }
                    break
                except subprocess.CalledProcessError as e:
                    logger.error("sacct process error")
                    logger.error(e)
                except IndexError as e:
                    pass
                # Try getting job with scontrol instead in case sacct is misconfigured
                if not res:
                    try:
                        sctrl_res = subprocess.check_output(
                            "scontrol -o show job {}".format(jobid), shell=True
                        )
                        m = re.search("JobState=(\w+)", sctrl_res.decode())
                        res = {jobid: m.group(1)}
                        break
                    except subprocess.CalledProcessError as e:
                        logger.error("scontrol process error")
                        logger.error(e)

                if i >= STATUS_ATTEMPTS - 1:
                    raise WorkflowError("Unable to query job status for 10 times")

        return res[jobid]  # == status

    def _wait_for_jobs(self):
        # busy wait on job completion
        # This is only needed if your backend does not allow to use callbacks
        # for obtaining job status.
        fail_stati = (
            "BOOT_FAIL",
            "OUT_OF_MEMORY",
            "CANCELLED",
            "FAILED",
            "NODE_FAIL",
            "DEADLINE",
            "PREEMPTED",
            "TIMEOUT",
            "ERROR",
        )
        while True:
            # always use self.lock to avoid race conditions
            with self.lock:
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for j in active_jobs:
                status = self.job_status(j)
                if status == "COMPLETED":
                    j.callback(j.job)
                elif status in fail_stati:
                    self.print_job_error(j.job, msg="failed with status {status}")
                    j.error_callback(j.job)
                else:  # still running?
                    still_running.append(j)

            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()
