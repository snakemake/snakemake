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

    #default_jobscript = "slurm/slurm_jobstep.py"

    def __init__(
        self,
        workflow,
        dag,
        cores,
        jobname="snakejob_{name}_{jobid}", # name and SLURM_JOB_ID will be appended
        printreason=False,
        quiet=False,
        printshellcmds=False,
        latency_wait=3,
        regions=None,
        location=None,
        cache=False,
        local_input=None,
        restart_times=0,
        exec_job=None,
        max_status_checks_per_second=1,
        cluster_config=None,
    ):
        #self.jobname = self.dag.rulename + '_' + self.dag.jobid
        # needs to be set in either case for the submission string
        if not cores:
            cores = 1

        # seems otherswise not in namespace, despite the super class init
        self.max_status_checks_per_second = max_status_checks_per_second

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
            assume_shared_fs=True,
            max_status_checks_per_second=max_status_checks_per_second,
        )
        self.exec_job += " --slurm-jobstep"

    def cancel(self):
        for job in self.active_jobs:
            jobid = job.jobid
            try:
                # timeout set to 60, because a scheduler cycle usually is
                # about 30 sec, but can be longer in extreme cases.
                # Under 'normal' circumstances, 'scancel' is executed in
                # virtually no time.
                subprocess.run(shlex.split("scancel {}".format(jobid)), timeout=60)
            except subprocess.TimeoutExpired:
                pass  # shell we ignore the timeout, here?
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

    def cluster_params(self, job):
        """
           Returns wildcards object for 'job'.

           In contrast to the ClusterExecutor, which gets
           its config from a config file, this SlurmExecutor
           has the internal handling via job.resources
        """
        #print('in cluster_params():')
        #print('========')
        #print('type(job):')
        #print(type(job))
        #print('-------')
        #print('dir(job):')
        #print(dir(job))
        #print('-------')
        #cluster = job.dynamic_wildcards.copy()
        #print('in cluster_params():')
        ###print(cluster)
        #sys.exit()

        return job.dynamic_wildcards.copy()

    def write_jobscript(self, job, jobscript, **kwargs):
        pass
        # only force threads if this is not a group job
        # otherwise we want proper process handling
    #    use_threads = "--force-use-threads" if not job.is_group() else ""

    #    envvars = " ".join(
    #       "{}={}".format(var, os.environ[var]) for var in self.workflow.envvars
    #    )

        # We need to shorten the exec_job-string, here.
        # The first part SLURM sees should be the executable, here:
        # snakemake.
        #print(self.exec_job)
        #print('=' * 10)
        #to_exec = " ".join(self.exec_job.split(" ")[5:])
        #print(to_exec)
        #print('=' * 10)
        #print(kwargs)
        #print('=' * 10)
        #print(job)
        #print('=' * 10)
        #print(self.jobscript)
        #print('=' * 10)

    #    exec_job = self.format_job(
    #        #to_exec,
    #        self.exec_job,
    #        job,
    #        _quote_all=True,
    #        use_threads=use_threads,
    #        #envvars=envvars,
    #        **kwargs,
    #    )

        #content = self.format_job(self.jobscript, job, exec_job=exec_job, **kwargs)#

        #logger.debug("Jobscript:\n{}".format(content))
        #with open(jobscript, "w") as f:
        #    print(content, file=f)

        #os.chmod(jobscript, os.stat(jobscript).st_mode | stat.S_IXUSR | stat.S_IRUSR)

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        super()._run(job)
        workdir = os.getcwd()
        jobid = job.jobid
        # generic part of a submission string:
        # print(self.exec_job)
        os.makedirs('.snakemake/slurm_logs', exist_ok=True)
        #print('got here: run() in slurm_submit')
        #sys.exit()

        envvars = " ".join(
            "{}={}".format(var, os.environ[var]) for var in self.workflow.envvars
        )
        jobname = self.jobname.format()
        try:
            call = "sbatch -A {account} -p {partition} \
                    -J {jobname} \
                    -o .snakemake/slurm_logs/%x_%j.log \
                    --export=ALL".format(
                **job.resources, jobname=jobname
            )
        except KeyError as e:
            logger.error(
                "Missing job submission key '{}' for job '{}'.".format(
                    e.args[0], job.name
                )
            )
            sys.exit(1)
        if not job.resources.get('walltime_minutes'):
            logger.warning("No wall time limit is set, setting 'walltime_minutes' to 1.")
        call += " -t {walltime_minutes}".format(walltime_minutes = job.resources.get('walltime_minutes', default_value=1))
        if job.resources.get("constraint"):
            call += " -C {constraint}".format(**job.resources)
        # TODO: implement when tempfs-resource is defined
        # if job.resources.get()
        if job.resources.get("mem_mb_per_cpu"):
            call += " --mem-per-cpu={mem_mb_per_cpu}".format(**job.resources)
        elif job.resources.get("mem_mb"):
            call += " --mem {mem_mb}".format(**job.resources)
        else:
            logger.warning(
                "No job memory information ('mem_mb' or 'mem_mb_per_cpu') is given - submitting without. This might or might not work on your cluster."
            )

        # if ressources == MPI:
        if job.resources.get("mpi", False):
            if job.resources.get("nodes", False):
                call += " --nodes={}".format(job.resources.get("nodes", 1))
            if job.resources.get("tasks", False):
                call += " --ntasks={}".format(job.resources.get("tasks", 1))
            if job.resources.get("threads", False):
                call += " -c={}".format(job.resources.get("threads", 1))
        # ordinary smp application:
        elif not job.is_group():
            # TODO: this line will become longer
            # TODO: hence the single command, yet
            if job.threads == 0:
                call += " -n 1 -c 1 {exec_job}".format(self.exec_job)
            else:
                call += " -n 1 -c {threads}".format(threads=job.threads)
        # job.is_group
        else:
            ntasks = max(map(len, job.toposorted))
            threads = max(j.threads for j in job)
            call += " -n {ntasks} -c {threads}".format(ntasks=ntasks,
                                     threads=threads)


        # as we cannot do 'sbatch ... cd {workflow.workdir_init} && python ...' as snakemake otherwise expects this
        # usually, this is default SLURM behaviour, but to ensure this working
        # on all clusters, we make it explicit
        call += " --chdir={workdir_init}".format(
            workdir_init=self.workflow.workdir_init
        )
        # and finally the job to execute with all the snakemake parameters

        #jobscript = self.get_jobscript(job)
        jobfinished = os.path.join(self.tmpdir, "{}.jobfinished".format(jobid))
        jobfailed = os.path.join(self.tmpdir, "{}.jobfailed".format(jobid))
        #self.write_jobscript(
        #    job, jobscript, jobfinished=jobfinished, jobfailed=jobfailed
        #)

        exec_job = self.format_job(
            self.exec_job,
            job,
            envvars=envvars,
            _quote_all=True,
            use_threads="--force-use-threads" if not job.is_group() else "",
        )

        # we need to modify `exec_job`, because we can only wrap executables
        # and not `cd` statements or tinkerings with PATH:
        # the directory is changed anyways by SLURM, the environment is
        # carried on, within envvars and the SLURM environment.
        #exe
        #exec_job = ' '.join(exec_job.split(' ')[8:]).replace('\\', ' ').replace(os.linesep, ' ')
        index = exec_job.index("-m snakemake")
        exec_job=exec_job[index:]
        call += ' --wrap=\'python {exec_job} --jobs unlimited\''.format(exec_job=exec_job)
        #sys.exit()
        #try:
        out = subprocess.check_output(call, shell=True, encoding="ascii").strip()
        #except:
        #    pass  # check template

        jobid = out.split(" ")[-1]
        self.active_jobs.append(SlurmJob(job, jobid, callback, error_callback))

    def job_status(self, jobid: int):
        """
        obtain SLURM job status of submitted jobs
        """
        STATUS_ATTEMPS = 10
        res = None

        for i in range(STATUS_ATTEMPS):
            # use self.status_rate_limiter to avoid too many API calls.
            with self.status_rate_limiter:
                try:
                    sacct_cmd = shlex.split("sacct -P -b -j {} -n".format(jobid))
                    sacct_res = subprocess.check_output(sacct_cmd, encoding="ascii")
                    res = {
                        x.split("|")[0]: x.split("|")[1]
                        for x in sacct_res.strip().split("\n")
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
                        sctrl_cmd = shlex.split("scontrol -o show job {}".format(jobid))
                        sctrl_res = subprocess.check_output(sctrl_cmd, encoding="ascii")
                        m = re.search("JobState=(\w+)", sctrl_res)
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
                status = self.job_status(j.jobid)
                if status == "COMPLETED":
                    j.callback(j.job)
                elif status in fail_stati:
                    self.print_job_error(j.job, msg="failed with SLURM status '{status}'".format(status=status))
                    j.error_callback(j.job)
                else:  # still running?
                    still_running.append(j)

            with self.lock:
                self.active_jobs.extend(still_running)
            time.sleep(1 / self.max_status_checks_per_second)
