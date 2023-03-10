from collections import namedtuple
from functools import partial
from io import StringIO
import os
import re
import stat
import sys
import time
import shlex
import shutil
import subprocess
import tarfile
import tempfile
import uuid
import pandas as pd

from snakemake.jobs import Job
from snakemake.logging import logger
from snakemake.exceptions import print_exception
from snakemake.exceptions import log_verbose_traceback
from snakemake.exceptions import WorkflowError
from snakemake.executors import ClusterExecutor
from snakemake.utils import makedirs
from snakemake.io import get_wildcard_names, Wildcards
from snakemake.common import async_lock

SlurmJob = namedtuple("SlurmJob", "job jobid callback error_callback slurm_logfile")


def get_account():
    """
    tries to deduce the acccount from recent jobs,
    returns None, if none is found
    """
    cmd = f'sacct -nu "{os.environ["USER"]}" -o Account%20 | head -n1'
    try:
        sacct_out = subprocess.check_output(
            cmd, shell=True, text=True, stderr=subprocess.PIPE
        )
        return sacct_out.strip()
    except subprocess.CalledProcessError as e:
        logger.warning(
            f"No account was given, not able to get a SLURM account via sacct: {e.stderr}"
        )
        return None


def test_account(account):
    """
    tests whether the given account is registered, raises an error, if not
    """
    cmd = f'sacctmgr -n -s list user "{os.environ["USER"]}" format=account%20'
    try:
        accounts = subprocess.check_output(
            cmd, shell=True, text=True, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise WorkflowError(
            f"Unable to test the validity of the given or guessed SLURM account '{account}' with sacctmgr: {e.stderr}"
        )

    accounts = accounts.split()

    if account not in accounts:
        raise WorkflowError(
            f"The given account {account} appears to be invalid. Available accounts:\n{', '.join(accounts)}"
        )


def get_default_partition(job):
    """
    if no partition is given, checks whether a fallback onto a default partition is possible
    """
    try:
        out = subprocess.check_output(
            r"sinfo -o %P", shell=True, text=True, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise WorkflowError(
            f"Failed to run sinfo for retrieval of cluster partitions: {e.stderr}"
        )
    for partition in out.split():
        # a default partition is marked with an asterisk, but this is not part of the name
        if "*" in partition:
            # the decode-call is necessary, because the output of sinfo is bytes
            return partition.replace("*", "")
    logger.warning(
        f"No partition was given for rule '{job}', and unable to find a default partition."
        " Trying to submit without partition information."
        " You may want to invoke snakemake with --default-resources 'slurm_partition=<your default partition>'."
    )
    return ""


class SlurmExecutor(ClusterExecutor):
    """
    the SLURM_Executor abstracts execution on SLURM
    clusters using snakemake resource string
    """

    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob_{name}_{jobid}",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        restart_times=0,
        max_status_checks_per_second=2,
        cluster_config=None,
    ):
        super().__init__(
            workflow,
            dag,
            cores,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            cluster_config=cluster_config,
            restart_times=restart_times,
            assume_shared_fs=True,
            max_status_checks_per_second=max_status_checks_per_second,
        )
        self.run_uuid=str(uuid.uuid4())
        self._fallback_account_arg = None
        self._fallback_partition = None

    def additional_general_args(self):
        # we need to set -j to 1 here, because the behaviour
        # of snakemake is to submit all jobs at once, otherwise.
        # However, the SLURM Executor is supposed to submit jobs
        # one after another, so we need to set -j to 1 for the
        # JobStep Executor, which in turn handles the launch of
        # SLURM jobsteps.
        return [" --slurm-jobstep", "--jobs 1"]

    def cancel(self):
        # Jobs are collected to reduce load on slurmctld
        jobids = " ".join([job.jobid for job in self.active_jobs])
        if len(jobids) > 0:
            try:
                # timeout set to 60, because a scheduler cycle usually is
                # about 30 sec, but can be longer in extreme cases.
                # Under 'normal' circumstances, 'scancel' is executed in
                # virtually no time.
                subprocess.check_output(
                    f"scancel {jobids}",
                    text=True,
                    shell=True,
                    timeout=60,
                    stderr=subprocess.PIPE,
                )
            except subprocess.TimeoutExpired:
                logger.warning(f"Unable to cancel jobs within a minute.")
        self.shutdown()

    def get_account_arg(self, job):
        """
        checks whether the desired account is valid,
        returns a default account, if applicable
        else raises an error - implicetly.
        """
        if job.resources.get("slurm_account"):
            # here, we check whether the given or guessed account is valid
            # if not, a WorkflowError is raised
            test_account(job.resources.slurm_account)
            return f" -A {job.resources.slurm_account}"
        else:
            if self._fallback_account_arg is None:
                logger.warning("No SLURM account given, trying to guess.")
                account = get_account()
                if account:
                    logger.warning(f"Guessed SLURM account: {account}")
                    self._fallback_account_arg = f" -A {account}"
                else:
                    logger.warning(
                        "Unable to guess SLURM account. Trying to proceed without."
                    )
                    self._fallback_account_arg = (
                        ""  # no account specific args for sbatch
                    )
            return self._fallback_account_arg

    def get_partition_arg(self, job):
        """
        checks whether the desired partition is valid,
        returns a default partition, if applicable
        else raises an error - implicetly.
        """
        if job.resources.get("slurm_partition"):
            partition = job.resources.slurm_partition
        else:
            if self._fallback_partition is None:
                self._fallback_partition = get_default_partition(job)
            partition = self._fallback_partition
        if partition:
            return f" -p {partition}"
        else:
            return ""

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        jobid = job.jobid

        log_folder = f"group_{job.name}" if job.is_group() else f"rule_{job.name}"

        slurm_logfile = f".snakemake/slurm_logs/{log_folder}/%j.log"
        os.makedirs(os.path.dirname(slurm_logfile), exist_ok=True)

        # generic part of a submission string:
        # we use a run_uuid as the job-name, to allow `--name`-based
        # filtering in the job status checks (`sacct --name` and `squeue --name`)
        call = f"sbatch --job-name {self.run_uuid} -o {slurm_logfile} --export=ALL"

        call += self.get_account_arg(job)
        call += self.get_partition_arg(job)

        if job.resources.get("runtime"):
            call += f" -t {job.resources.runtime}"
        else:
            logger.warning(
                "No wall time information given. This might or might not work on your cluster. "
                "If not, specify the resource runtime in your rule or as a reasonable "
                "default via --default-resources."
            )

        if job.resources.get("constraint"):
            call += f" -C {job.resources.constraint}"
        if job.resources.get("mem_mb_per_cpu"):
            call += f" --mem-per-cpu {job.resources.mem_mb_per_cpu}"
        elif job.resources.get("mem_mb"):
            call += f" --mem {job.resources.mem_mb}"
        else:
            logger.warning(
                "No job memory information ('mem_mb' or 'mem_mb_per_cpu') is given "
                "- submitting without. This might or might not work on your cluster."
            )

        # MPI job
        if job.resources.get("mpi", False):
            if job.resources.get("nodes", False):
                call += f" --nodes={job.resources.get('nodes', 1)}"
            if job.resources.get("tasks", False):
                call += f" --ntasks={job.resources.get('tasks', 1)}"

        cpus_per_task = job.threads
        if job.resources.get("cpus_per_task"):
            if not isinstance(cpus_per_task, int):
                raise WorkflowError(
                    "cpus_per_task must be an integer, but is {}".format(cpus_per_task)
                )
            cpus_per_task = job.resources.cpus_per_task
        # ensure that at least 1 cpu is requested
        # because 0 is not allowed by slurm
        cpus_per_task = max(1, cpus_per_task)
        call += f" --cpus-per-task={cpus_per_task}"

        if job.resources.get("slurm_extra"):
            call += f" {job.resources.slurm_extra}"

        exec_job = self.format_job_exec(job)
        # ensure that workdir is set correctly
        # use short argument as this is the same in all slurm versions
        # (see https://github.com/snakemake/snakemake/issues/2014)
        call += f" -D {self.workflow.workdir_init}"
        # and finally the job to execute with all the snakemake parameters
        call += f" --wrap={shlex.quote(exec_job)}"

        logger.debug(f"sbatch call: {call}")
        try:
            out = subprocess.check_output(
                call, shell=True, text=True, stderr=subprocess.STDOUT
            ).strip()
        except subprocess.CalledProcessError as e:
            raise WorkflowError(
                f"SLURM job submission failed. The error message was {e.output}"
            )

        slurm_jobid = out.split(" ")[-1]
        slurm_logfile = slurm_logfile.replace("%j", slurm_jobid)
        logger.info(
            f"Job {jobid} has been submitted with SLURM jobid {slurm_jobid} (log: {slurm_logfile})."
        )
        self.active_jobs.append(
            SlurmJob(job, slurm_jobid, callback, error_callback, slurm_logfile)
        )

    async def job_stati(self, command):
        """obtain SLURM job status of all submitted jobs with sacct

        Keyword arguments:
        command -- a slurm command that returns one line for each job with: 
                   "<raw/main_job_id>|<long_status_string>"
        """
        error = None
        try:
            time_before_query = time.time()
            command_res = subprocess.check_output(
                command, text=True, shell=True, stderr=subprocess.PIPE
            )
            query_duration = time.time() - time_before_query
            logger.debug(
                f"The job status was queried with command: {command}\n"
                f"It took: {query_duration} seconds\n"
                f"The output is:\n'{command_res}'\n"
            )
            res = dict(
                pd.read_csv(
                    StringIO(
                        command_res
                    )
                )
            )
        except subprocess.CalledProcessError as e:
            error = e.stderr
            pass
        except IndexError as e:
            pass

        return (res, error, query_duration)

    async def _wait_for_jobs(self):
        # busy wait on job completion
        # This is only needed if your backend does not allow to use callbacks
        # for obtaining job status.
        fail_stati = (
            "BOOT_FAIL",
            "CANCELLED",
            "DEADLINE",
            "FAILED",
            "NODE_FAIL",
            "OUT_OF_MEMORY",
            "PREEMPTED",
            "TIMEOUT",
            "ERROR",
        )
        # intialize time to sleep in seconds
        MIN_SLEEP_TIME = 20
        sleepy_time = MIN_SLEEP_TIME
        while True:
            # Initialize all query durations to specified 
            # 5 times the status_rate_limiter, to hit exactly
            # the status_rate_limiter for the first async below.
            # It is dynamically updated afterwards.
            sacct_query_duration = squeue_query_duration = (1/self.status_rate_limiter)*5
            # always use self.lock to avoid race conditions
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                active_jobs_ids = { j.jobid for j in active_jobs }
                self.active_jobs = list()
                still_running = list()
            STATUS_ATTEMPTS = 10
            # this code is inspired by the snakemake profile: TODO: link to github
            for i in range(STATUS_ATTEMPTS):
                # use self.status_rate_limiter and adaptive query
                # timing to avoid too many API calls.
                async with min(
                    self.status_rate_limiter,
                    # if slurmdbd (sacct) is strained and slow, reduce the query frequency
                    (1/sacct_query_duration)/5,
                    # if slurmctld (squeue) is strained and slow, reduce the query frequency
                    (1/squeue_query_duration)/5,
                ):
                    (status_of_jobs, sacct_err, sacct_query_duration) = await self.job_stati(
                        # -X: only show main job, no substeps
                        f"sacct -X --parsable2 --noheader --format=JobIdRaw,State --name {self.run_uuid}"
                    )
                    no_sacct_status = active_jobs_ids - set(status_of_jobs.keys())
                    if not no_sacct_status:
                        break
                    else:
                        job_ids_string = ",".join(no_sacct_status)
                        (squeue_stati, squeue_err, squeue_query_duration) = await self.job_stati(
                            # * %F: job array ID gives either the main ID for the job array,
                            #       or the regular job id for non-array jobs (this should be
                            #       equivalent to JobIdRaw in sacct)
                            # * %T: long job status
                            f"squeue --noheader --format=\"%F|%T\" --name {self.run_uuid} --jobs {job_ids_string}"
                        )
                        status_of_jobs.update(squeue_stati)

                        remaining_without_status = active_jobs_ids - status_of_jobs
                        if not remaining_without_status:
                            break
                        else:
                            def fmt_err(err_type, err_msg):
                                if err_msg is not None:
                                    return f"  {err_type} error: {err_msg.strip()}"
                                else:
                                    return ""

                            logger.error(
                                f"Could not get status of slurm jobs:\n"
                                f"  {','.join(remaining_without_status)}\n"
                                f"If the respective status queries reported errors, those were:\n"
                                f"{fmt_err('sacct', sacct_err)}\n"
                                f"{fmt_err('scontrol', squeue_err)}\n"
                            )
                if i >= STATUS_ATTEMPTS - 1:
                    raise WorkflowError(f"Unable to get the status of all active_jobs with sacct and squeue, even after {STATUS_ATTEMPTS} attempts")
            for j in active_jobs:
                status = status_of_jobs[j.jobid]
                if status == "COMPLETED":
                    j.callback(j.job)
                elif status == "UNKNOWN":
                    # the job probably does not exist anymore, but 'sacct' did not work
                    # so we assume it is finished
                    j.callback(j.job)
                elif status in fail_stati:
                    self.print_job_error(
                        j.job,
                        msg=f"SLURM-job '{j.jobid}' failed, SLURM status is: '{status}'",
                        aux_logs=[j.slurm_logfile],
                    )
                    j.error_callback(j.job)
                else:  # still running?
                    still_running.append(j)

            # no jobs finished in the last query period
            if not active_jobs_ids - { j.jobid for j in still_running }:
                # sleep a little longer, but never longer than a max
                sleepy_time = min(sleepy_time + 20, 240)
            else:
                sleepy_time = MIN_SLEEP_TIME
            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            time.sleep(sleepy_time)
