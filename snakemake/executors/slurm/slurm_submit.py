from collections import namedtuple
from functools import partial
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
        max_status_checks_per_second=0.03,
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
        for job in self.active_jobs:
            jobid = job.jobid
            try:
                # timeout set to 60, because a scheduler cycle usually is
                # about 30 sec, but can be longer in extreme cases.
                # Under 'normal' circumstances, 'scancel' is executed in
                # virtually no time.
                subprocess.check_output(
                    f"scancel {jobid}",
                    text=True,
                    shell=True,
                    timeout=60,
                    stderr=subprocess.PIPE,
                )
            except subprocess.TimeoutExpired:
                logger.warning(f"Unable to cancel job {jobid} within a minute.")
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
        call = f"sbatch -J {self.get_jobname(job)} -o {slurm_logfile} --export=ALL"

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

    async def job_status(self, jobid: int):
        """
        obtain SLURM job status of submitted jobs
        """
        STATUS_ATTEMPTS = 10
        res = None
        # this code is inspired by the snakemake profile: TODO: link to github
        for i in range(STATUS_ATTEMPTS):
            # use self.status_rate_limiter to avoid too many API calls.
            async with self.status_rate_limiter:
                sacct_error = None
                try:
                    sacct_cmd = f"sacct -P -n --format=JobIdRaw,State -j {jobid}"
                    sacct_res = subprocess.check_output(
                        sacct_cmd, text=True, shell=True, stderr=subprocess.PIPE
                    )
                    logger.debug(f"The sacct output is: '{sacct_res}'")
                    res = {
                        x.split("|")[0]: x.split("|")[1]
                        for x in sacct_res.strip().split("\n")
                    }
                    break
                except subprocess.CalledProcessError as e:
                    sacct_error = e.stderr
                    pass  # try scontrol below
                except IndexError as e:
                    pass
                # Try getting job with scontrol instead in case sacct is misconfigured
                if not res:
                    try:
                        sctrl_cmd = f"scontrol show jobid -dd {jobid}"
                        out = subprocess.check_output(
                            sctrl_cmd,
                            shell=True,
                            stderr=subprocess.PIPE,
                            text=True,
                        )
                        logger.debug(f"The scontrol output is: '{out}'")
                        m = re.search(r"JobState=(\w+)", out)
                        res = {jobid: m.group(1)}
                        break
                    except subprocess.CalledProcessError as e:

                        def fmt_err(err_type, err_msg):
                            if err_msg is not None:
                                return f"\n    {err_type} error: {err_msg.strip()}"
                            else:
                                return ""

                        logger.error(
                            f"Error getting status of slurm job {jobid}:{fmt_err('sacct', sacct_error)}{fmt_err('scontrol', e.stderr)}"
                        )

                if i >= STATUS_ATTEMPTS - 1:
                    raise WorkflowError("Unable to query job status for 10 times")

        return res[jobid]  # == status

    async def _wait_for_jobs(self):
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
            async with async_lock(self.lock):
                if not self.wait:
                    return
                active_jobs = self.active_jobs
                self.active_jobs = list()
                still_running = list()
            for j in active_jobs:
                status = await self.job_status(j.jobid)
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

            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)
            time.sleep(1 / self.max_status_checks_per_second)
