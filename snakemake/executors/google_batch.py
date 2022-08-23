from collections import namedtuple


GoogleLifeSciencesJob = namedtuple(
    "GoogleLifeSciencesJob", "job jobname jobid callback error_callback"
)


class SkeletonExecutor(ClusterExecutor):
    def __init__(
        self,
        workflow,
        dag,
        jobname="snakejob.{name}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        latency_wait=3,
        cluster_config=None,
        local_input=None,
        restart_times=None,
        max_status_checks_per_second=1,
    ):

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
            assume_shared_fs=False,
            max_status_checks_per_second=max_status_checks_per_second,
        )

        # add additional attributes

    def shutdown(self):
        # perform additional steps on shutdown if necessary
        super().shutdown()

    def cancel(self):
        for job in self.active_jobs:
            # cancel active jobs here
            pass
        self.shutdown()

    def run_jobs(self, jobs, callback=None, submit_callback=None, error_callback=None):
        """Run a list of jobs that is ready at a given point in time.

        By default, this method just runs each job individually.
        This behavior is inherited and therefore this method can be removed from the skeleton if the
        default behavior is intended.
        This method can be overwritten to submit many jobs in a more efficient way than one-by-one.

        Note that in any case, for each job, the callback functions have to be called individually!
        """
        for job in jobs:
            self.run(
                job,
                callback=callback,
                submit_callback=submit_callback,
                error_callback=error_callback,
            )

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        """Run an individual job or a job group."""

        # Necessary: perform additional executor independent steps before running the job
        super()._run(job)

        # obtain job execution command
        exec_job = self.format_job(
            self.exec_job,
            job,
            _quote_all=True,
            use_threads="--force-use-threads" if not job.is_group() else "",
        )

        # submit job here, and obtain job ids from the backend

        # register job as active, using your own namedtuple.
        # The namedtuple must at least contain the attributes
        # job, jobid, callback, error_callback.
        self.active_jobs.append(
            GoogleLifeSciencesJob(job, jobid, callback, error_callback)
        )

    def _wait_for_jobs(self):
        from snakemake.executors import sleep

        # busy wait on job completion
        # if possible, use batch requests for many jobs, to avoid problems
        # with API rate limits.
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

                    # Retrieve status of job j from your backend via j.jobid
                    # Handle completion and errors, calling either j.callback(j.job)
                    # or j.error_callback(j.job)
                    # In case of error, add job j to still_running.
                    pass
            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()
