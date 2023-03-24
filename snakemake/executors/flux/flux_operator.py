__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import json
import os
import re
import time
from collections import namedtuple

from snakemake.executors import ClusterExecutor, sleep
from snakemake.logging import logger
from snakemake.resources import DefaultResources
from snakemake.common import async_lock
from snakemake.exceptions import WorkflowError

# Import Flux Operator SDK
try:
    import fluxoperator
    from fluxoperator.client import FluxMiniCluster
except ImportError:
    fluxoperator = None
    FluxMiniCluster = object

# uid is for the MiniCluster
FluxJob = namedtuple(
    "FluxJob", "job jobname jobid callback error_callback uid is_checkpoint"
)
JobInfo = namedtuple("JobInfo", "jobid user name state ntasks nnodes time node")

# completed, failed, cancelled, timeout
done_states = ["CD", "F", "CA", "TO"]

# Shared flux user
fluxuser = "fluxuser"
socket = "local:///run/flux/local"

# Wrap the MiniCluster to add extra execution logic
class MiniClusterExec(FluxMiniCluster):

    def execute(self, command, print_result=False, quiet=True):
        """
        Wrap the kubectl_exec to add logic to issue to the broker instance.
        """
        envars = ["PATH", "PYTHONPATH", "LD_LIBRARY_PATH"]
        envars = " ".join([f"-E {e}=${e}" for e in envars])
        home = f"/home/{fluxuser}"

        # Always execute to the broker pod
        res = self.ctrl.kubectl_exec(
            f"sudo -u {fluxuser} {envars} -E HOME={home} flux proxy {socket} {command}",
            quiet=quiet,
            pod=self.broker_pod,
            namespace=self.namespace,
            name=self.name,
        )
        if print_result:
            print(res, end="")
        return res

    def wait_for_broker_socket(self):
        """
        Wrap the kubectl_exec to add logic to issue to the broker instance.
        """
        while True:
            res = self.execute("flux jobs -a")
            logger.debug(res)
            # flux-jobs: ERROR: No service matching job-list.list is registered
            if "JOBID" in res:
                break
            time.sleep(3)


class FluxManager:
    """
    The manager keeps track of Flux MiniClusters.
    """

    def __init__(self, namespace, *args, **kwargs):
        self.namespace = namespace

        # Keep track of active MiniClusters, assignments
        self._clusters = {}
        self.assignment = {}

        # This will be the job that (after completion) we can delete a minicluster
        # This is essentially when a uid is last seen
        self._deletions = {}

    def check_deletion(self, jobname, uid):
        """
        When a job finishes, check if the MiniCluster is good to delete.
        """
        if self._deletions.get(uid) == jobname:
            self.cleanup_miniclusters(uid)

    def jobinfo(self, jobid, uid):
        """
        Get job info for a jobid
        """
        raw = self._clusters[uid].execute(f"flux jobs {jobid} --no-header")
        parts = [x for x in re.split("[ ]+", raw) if x]
        return JobInfo(*parts)

    def job_log(self, jobid, uid):
        """
        Get the log for a failed job (this will hang otherwise)
        """
        return self._clusters[uid].execute(f"flux job attach {jobid}")

    def has_cluster(self, name):
        """
        Return a MiniCluster by name.
        """
        return self.get_cluster(name) is not None

    def reset_assignments(self):
        """
        Reset assignments of MiniClusters to jobs.
        """
        self.assignment = {}
        self._deletions = {}

    def set_deletion_after(self, uid, jobname):
        """
        Keep track of when a MiniCluster can be deleted.

        This is the last job that needs it.
        """
        self._deletions[uid] = jobname

    def assign(self, jobname, nodes, image):
        """
        Assign a job name to a MiniCluster

        This uses the nodes and image to generate a unique id.
        """
        slug = re.sub("([]:]|[/]|[.])", "-", image)
        uid = f"{slug}-{nodes}"
        self.assignment[jobname] = {"nodes": nodes, "uid": uid, "image": image}
        return uid

    def get_cluster(self, name):
        """
        Return a MiniCluster by name.
        """
        return self._clusters.get(name)

    def cancel_job(self, jobid):
        """
        Cancel a running job
        """
        res = self.execute("flux job cancel {jobid}")
        logger.debug(res)

    def cleanup_miniclusters(self, name=None):
        """
        Delete running MiniClusters
        """
        if name and name not in self._clusters:
            logger.warning(f"WARNING: {name} is not a known MiniCluster")
            return

        import IPython
        IPython.embed()
        to_delete = [name] if name else list(self._clusters.keys())
        for name in to_delete:
            logger.info(f"Deleting MiniCluster {name}")
            self._clusters[name].delete()
            del self._clusters[name]

    def new_minicluster(self, minicluster, container):
        """
        Create a new MiniCluster and wait for it to be ready.
        """
        logger.info(f"Creating new MiniCluster {minicluster['name']}...")

        # Add the shared namespace
        minicluster["namespace"] = self.namespace

        # The operator will time creation through pods being ready
        mc = MiniClusterExec()
        mc.create(**minicluster, container=container)

        # Add a new MiniCluster - the name is the uid we had previously
        self._clusters[mc.name] = mc

        # Give time for broker to start (this likely needs to vary)
        logger.debug("Giving broker time to start...")
        mc.wait_for_broker_socket()
        return mc


class FluxOperatorExecutor(ClusterExecutor):
    """
    The Flux operator executor deploys jobs each to a Flux MiniCluster on Kubernetes

    The max nodes (the size of the MiniCluster) could be set via the client init, but
    for now is expected to be defined in the environment (or default to 4 for tests).
    """

    def __init__(
        self,
        workflow,
        dag,
        jobname="snakejob.{name}.{jobid}.sh",
        printreason=False,
        quiet=False,
        printshellcmds=False,
        kubernetes_nodes=4,
        flux_operator_ns="flux-operator",
        container_image=None,
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

        # Set the default container image
        self.container_image = container_image or 'ghcr.io/rse-ops/mamba:app-mamba'

        # Attach variables for easy access
        self.workdir = os.path.realpath(os.path.dirname(self.workflow.persistence.path))
        self.envvars = list(self.workflow.envvars) or []
        self._plan = {}

        # Quit early if we can't access the flux api
        if not fluxoperator:
            raise WorkflowError(
                "Please install the fluxoperator python package (pip install fluxoperator)"
            )
        self.ctrl = FluxManager(flux_operator_ns)

        # Set max nodes of k8s cluster
        self.set_max_nodes(kubernetes_nodes)

        # Plan execution of dag - minicluster assignments
        # This also will be run when a checkpoint is completed - job.is_checkpoint is True
        self.update_plan()

    def shutdown(self):
        """
        Shutdown and cleanup MiniClusters
        """
        super().shutdown()
        self.ctrl.cleanup_miniclusters()

    def set_max_nodes(self, max_nodes):
        """
        Set the max nodes for the MiniCluster (typically the k8s cluster size)
        """
        # Set a default of 8 max nodes, unless in the environment
        try:
            self._max_nodes = int(os.environ.get("SNAKEMAKE_FLUX_OPERATOR_MAX_NODES"))
        except:
            self._max_nodes = max_nodes

    def cancel(self):
        """
        cancel execution, usually by way of control+c. Cleanup is done in
        shutdown (deleting cached workdirs in Google Cloud Storage
        """
        for job in self.active_jobs:
            jobinfo = self.ctrl.jobinfo(job.jobid, job.uid)
            if jobinfo.state not in done_states:
                self.ctrl.cancel_job(job.jobid)
        self.shutdown()

    def update_plan(self):
        """
        Derive an assignment of each job level to a MiniCluster size.

        Assuming that the jobs can run together / at once, we can run
        them on the same MiniCluster, and can derive their needed size.
        This also takes into account the image. This assumes that jobs
        on the same level are similar in size. If they are not, we will
        need to check for some range of difference (and then assign
        to different MiniClusters).
        """
        max_nodes = self._max_nodes

        # Reset current assignments
        self.ctrl.reset_assignments()

        # Each job will be assigned a number of nodes based on its level
        for level in self.dag.toposorted():
            # This is nodes needed for an entire level of jobs
            needed = 0
            for job in level:
                job_nodes = job.resources.get("_nodes", 1)

                # This cannot be satisfied if the cluster isn't large enough for the job
                if job_nodes > max_nodes:
                    raise WorkflowError(
                        "Job {job} requires {job_nodes}, but max nodes set to {max_nodes}"
                    )
                needed += job.resources.get("_nodes", 1)

            # Never create a size 1 MiniCluster, only a broker
            if needed == 1:
                needed += 1

            # We don't want to blow up a MiniCluster size beyond what is reasonable.
            # Don't ever go above the max size.
            if needed > max_nodes:
                logger.warning(
                    f"Level {level} requires {needed}, over the limit of {max_nodes}"
                )
                max_nodes = needed

            # Make the assignment
            for job in level:
                # All isn't normally an actual step
                if job.name == "all":
                    continue

                uid = self.ctrl.assign(job.name, needed, self.container_image)

                # The deletion for the UID happens at the time the job is last seen
                # When this job completes
                self.ctrl.set_deletion_after(uid, job.name)

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

    def ensure_minicluster(self, job):
        """
        Ensure a MiniCluster is running for a specific job.
        """
        # Get the assignment for the job
        assignment = self.ctrl.assignment[job.name]
        uid = assignment["uid"]

        # We already created the MiniCluster - we're good
        if self.ctrl.has_cluster(uid):
            return self.ctrl.get_cluster(uid), uid

        # TODO If we don't have the MiniCluster, ensure we have enough room for it
        # If we do, we can create, otherwise we need to cleanup or wait
        # The container has a named volume "data" bound at the working directory path
        container = {
            "image": assignment["image"],
            "run_flux": True,
            "volumes": {"data": {"path": self.workdir}},
            "flux_user": {"name": fluxuser},
        }

        # The MiniCluster is expecting /tmp/workflow to be bound on the node
        minicluster = {
            "namespace": self.ctrl.namespace,
            "interactive": True,
            "volumes": {"data": {"path": "/tmp/workflow", "storage_class": "hostpath"}},
            "name": uid,
            "size": assignment["nodes"],
        }
        logger.info(
            f"To delete this cluster: kubectl delete -n {self.ctrl.namespace} pod {uid}"
        )
        return self.ctrl.new_minicluster(minicluster, container), uid

    def run(self, job, callback=None, submit_callback=None, error_callback=None):
        """
        Submit a job to the Flux Operator MiniCluster
        """
        mc, uid = self.ensure_minicluster(job)
        super()._run(job)

        # Prepare job resourcces
        self._set_job_resources(job)

        # The entire snakemake command to run, etc
        command = self.format_job_exec(job)
        logger.debug(command)

        # Launch the flux job directly to cluster
        cores = job.resources.get("_cores", 1)
        nodes = job.resources.get("_nodes", 1)
        runtime = job.resources.get("runtime", 0)

        # Note that we can't include the current host environment
        res = mc.execute(
            f"flux submit --job-name {job} -t {runtime} --cores {cores} --nodes {nodes} --cwd {self.workdir} {command}"
        )
        jobid = res.strip()

        # A successful submit returns a job id
        if "ƒ" not in jobid:
            raise WorkflowError(f"There was an error submitting {job}: {jobid}")

        logger.info(f"Job {job.jobid} has been submitted with flux jobid {jobid}")

        # Waiting for the jobid is a small performance penalty, same as calling flux.job.submit
        self.active_jobs.append(
            FluxJob(
                job,
                job.name,
                jobid,
                callback,
                error_callback,
                uid,
                job.is_checkpoint,
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
                try:
                    jobinfo = self.ctrl.jobinfo(j.jobid, j.uid)
                except Exception as e:
                    self.print_job_error(j.job, jobid=j.jobid)
                    raise WorkflowError(f"Issue getting job info: {e}")

                # completed, failed, cancelled, timeout
                if jobinfo.state in done_states:
                    # If we hit the last job for a given minicluster
                    self.ctrl.check_deletion(j.jobname, j.uid)

                    # the job finished (but possibly with nonzero exit code)
                    if jobinfo.state == "F":
                        log = self.ctrl.job_log(j.jobid, j.uid)
                        logger.error(log)
                        self.print_job_error(j.job, jobid=j.jobid)
                        j.error_callback(j.job)
                        continue

                    # If it's a checkpoint, update the minicluster plan
                    if j.is_checkpoint:
                        self.update_plan()

                    # Finished and success!
                    j.callback(j.job)

                # Otherwise, we are still running
                else:
                    still_running.append(j)
            async with async_lock(self.lock):
                self.active_jobs.extend(still_running)

            # Sleeps for 10 seconds
            await sleep()
