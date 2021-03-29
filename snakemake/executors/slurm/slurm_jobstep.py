import os
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


class SlurmJobstepExecutor(ClusterExecutor):
    """
    executes SLURM jobsteps and is *only* instaniated in
    a SLURM job context
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
        latency_wait=3,
        cluster_config=None,
        local_input=None,
        restart_times=None,
        exec_job=None,
        assume_shared_fs=True,
        max_status_checks_per_second=1,
    ):

        # overwrite the command to execute a single snakemake job if necessary
        # exec_job = "..."

        super().__init__(
            workflow,
            dag,
            None,
            jobname=jobname,
            printreason=printreason,
            quiet=quiet,
            printshellcmds=printshellcmds,
            latency_wait=latency_wait,
            cluster_config=cluster_config,
            local_input=local_input,
            restart_times=restart_times,
            exec_job=exec_job,
            assume_shared_fs=False,
            max_status_checks_per_second=10,
        )
        self.mem_per_node = os.getenv("SLURM_MEM_PER_NODE")
        self.cpus_on_node = os.getenv("SLURM_CPUS_ON_NODE")
        self.jobid = os.getenv("SLURM_JOB_ID")

    def format_job(self, job):
        return super().self.format_job(
            self.exec_job, job, _quote_all=True, use_threads=True
        )

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)

        jobsteps = dict()

        if job.is_group():

            def get_call(level_job, level_id, aux=""):
                # we need this calculation, because of srun's greediness and
                # SLURM's limits: it is not able to limit the memory if not per cpu
                mem_per_cpu = level_job.resources.mem_mb // level_job.threads
                return (
                    "srun -J {level_id} --jobid {jobid}"
                    "--mem-per-cpu {mem_per_cpu} -c {level_job.threads}"
                    "--exclusive -n 1 {aux} {exec_job}".format(
                        level_id=level_id,
                        level_job=level_job,
                        aux=aux,
                        jobid=self.jobid,
                        exec_job=self.format_job(level_job),
                    )
                )

            for level in job.toposorted:
                for level_id, level_job in enumerate(level[:-1]):
                    # TODO: spice with srun template
                    jobsteps[level_job] = subprocess.Popen(
                        get_call(level_job, level_id), shell=True
                    )
                # this way, we ensure that level jobs depending on the current level get started
                jobsteps[level_job] = subprocess.Popen(
                    get_call(level_job, level_id + 1, aux="--singleton"), shell=True
                )
        else:  # we combine smp and MPI jobs in one case
            # TODO: if MPI-job append '--cpu_bind=rank'
            call = "srun --cpu_bind=q {exec_job}".format(exec_job=self.format_job(job))
            jobsteps[job] = subprocess.Popen(call, shell=True)

        # wait until all steps are finished
        return_codes = list()
        error = False
        for job, proc in jobsteps.items():
            if proc.wait() != 0:
                self.print_job_error(job)
                error = True
        if error:
            error_callback(job)
        else:
            callback(job)
