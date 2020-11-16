__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2020, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import logging
import os
import sys
import time
import shutil
import tarfile
import tempfile
from collections import namedtuple
import uuid
import re
import math

from snakemake.logging import logger
from snakemake.exceptions import print_exception
from snakemake.exceptions import log_verbose_traceback
from snakemake.exceptions import WorkflowError
from snakemake.executors import ClusterExecutor, sleep
from snakemake.common import get_container_image, get_file_hash

# https://github.com/googleapis/google-api-python-client/issues/299#issuecomment-343255309
logging.getLogger("googleapiclient.discovery_cache").setLevel(logging.ERROR)

GoogleLifeSciencesJob = namedtuple(
    "GoogleLifeSciencesJob", "job jobname jobid callback error_callback"
)


class GoogleLifeSciencesExecutor(ClusterExecutor):
    """the GoogleLifeSciences executor uses Google Cloud Storage, and
    Compute Engine paired with the Google Life Sciences API.
    https://cloud.google.com/life-sciences/docs/quickstart
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
        container_image=None,
        regions=None,
        location=None,
        cache=False,
        latency_wait=3,
        local_input=None,
        restart_times=None,
        exec_job=None,
        max_status_checks_per_second=1,
        preemption_default=None,
        preemptible_rules=None,
    ):

        # Attach variables for easy access
        self.workflow = workflow
        self.quiet = quiet
        self.workdir = os.path.dirname(self.workflow.persistence.path)
        self._save_storage_cache = cache

        # Relative path for running on instance
        self._set_snakefile()

        # Prepare workflow sources for build package
        self._set_workflow_sources()

        exec_job = (
            exec_job
            or (
                "snakemake {target} --snakefile %s "
                "--force -j{cores} --keep-target-files --keep-remote "
                "--latency-wait {latency_wait} --scheduler {workflow.scheduler_type} "
                "--attempt 1 {use_threads} --max-inventory-time 0 "
                "{overwrite_config} {rules} --nocolor "
                "--notemp --no-hooks --nolock " % self.snakefile
            )
            + self.get_set_threads_args()
            + self.get_set_scatter_args()
        )

        # Set preemptible instances
        self._set_preemptible_rules(preemption_default, preemptible_rules)

        # IMPORTANT: using Compute Engine API and not k8s == no support secrets
        self.envvars = list(self.workflow.envvars) or []

        # Quit early if we can't authenticate
        self._get_services()
        self._get_bucket()

        # Akin to Kubernetes, create a run namespace, default container image
        self.run_namespace = str(uuid.uuid4())
        self.container_image = container_image or get_container_image()
        self.regions = regions or ["us-east1", "us-west1", "us-central1"]

        # The project name is required, either from client or environment
        self.project = (
            os.environ.get("GOOGLE_CLOUD_PROJECT") or self._bucket_service.project
        )

        # Determine API location based on user preference, and then regions
        self._set_location(location)

        # Tell the user right away the regions, location, and container
        logger.debug("regions=%s" % self.regions)
        logger.debug("location=%s" % self.location)
        logger.debug("container=%s" % self.container_image)

        # Keep track of build packages to clean up shutdown, and generate
        self._build_packages = set()
        targz = self._generate_build_source_package()
        self._upload_build_source_package(targz)

        # Save default resources to add later, since we need to add custom
        # default resources depending on the instance requested
        self.default_resources = self.workflow.default_resources
        self.workflow.default_resources.args = None

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

    def _get_services(self):
        """use the Google Discovery Build to generate API clients
        for Life Sciences, and use the google storage python client
        for storage.
        """
        from googleapiclient.discovery import build as discovery_build
        from oauth2client.client import (
            GoogleCredentials,
            ApplicationDefaultCredentialsError,
        )
        from google.cloud import storage

        # Credentials must be exported to environment
        try:
            creds = GoogleCredentials.get_application_default()
        except ApplicationDefaultCredentialsError as ex:
            log_verbose_traceback(ex)
            raise ex

        # Discovery clients for Google Cloud Storage and Life Sciences API
        self._storage_cli = discovery_build(
            "storage", "v1", credentials=creds, cache_discovery=False
        )
        self._compute_cli = discovery_build(
            "compute", "v1", credentials=creds, cache_discovery=False
        )
        self._api = discovery_build(
            "lifesciences", "v2beta", credentials=creds, cache_discovery=False
        )
        self._bucket_service = storage.Client()

    def _get_bucket(self):
        """get a connection to the storage bucket (self.bucket) and exit
        if the name is taken or otherwise invalid.

        Parameters
        ==========
        workflow: the workflow object to derive the prefix from
        """
        import google

        # Hold path to requested subdirectory and main bucket
        bucket_name = self.workflow.default_remote_prefix.split("/")[0]
        self.gs_subdir = re.sub(
            "^{}/".format(bucket_name), "", self.workflow.default_remote_prefix
        )
        self.gs_logs = os.path.join(self.gs_subdir, "google-lifesciences-logs")

        # Case 1: The bucket already exists
        try:
            self.bucket = self._bucket_service.get_bucket(bucket_name)

        # Case 2: The bucket needs to be created
        except google.cloud.exceptions.NotFound:
            self.bucket = self._bucket_service.create_bucket(bucket_name)

        # Case 2: The bucket name is already taken
        except Exception as ex:
            logger.error(
                "Cannot get or create {} (exit code {}):\n{}".format(
                    bucket_name, ex.returncode, ex.output.decode()
                )
            )
            log_verbose_traceback(ex)
            raise ex

        logger.debug("bucket=%s" % self.bucket.name)
        logger.debug("subdir=%s" % self.gs_subdir)
        logger.debug("logs=%s" % self.gs_logs)

    def _set_location(self, location=None):
        """The location is where the Google Life Sciences API is located.
        This can be meaningful if the requester has data residency
        requirements or multi-zone needs. To determine this value,
        we first use the locations API to determine locations available,
        and then compare them against:

        1. user specified location or prefix
        2. regions having the same prefix
        3. if cannot be satisifed, we throw an error.
        """
        # Derive available locations
        # See https://cloud.google.com/life-sciences/docs/concepts/locations
        locations = (
            self._api.projects()
            .locations()
            .list(name="projects/{}".format(self.project))
            .execute()
        )

        locations = {x["locationId"]: x["name"] for x in locations.get("locations", [])}

        # Alert the user about locations available
        logger.debug("locations-available:\n%s" % "\n".join(locations))

        # If no locations, there is something wrong
        if not locations:
            raise WorkflowError("No locations found for Google Life Sciences API.")

        # First pass, attempt to match the user-specified location (or prefix)
        if location:
            if location in locations:
                self.location = locations[location]
                return

            # It could be that a prefix was provided
            for contender in locations:
                if contender.startswith(location):
                    self.location = locations[contender]
                    return

            # If we get here and no match, alert user.
            raise WorkflowError(
                "Location or prefix requested %s is not available." % location
            )

        # If we get here, we need to select location from regions
        for region in self.regions:
            if region in locations:
                self.location = locations[region]
                return

        # If we get here, choose based on prefix
        prefixes = set([r.split("-")[0] for r in self.regions])
        regexp = "^(%s)" % "|".join(prefixes)
        for location in locations:
            if re.search(regexp, location):
                self.location = locations[location]
                return

        # If we get here, total failure of finding location
        raise WorkflowError(
            " No locations available for regions!"
            " Please specify a location with --google-lifesciences-location "
            " or extend --google-lifesciences-regions to find a Life Sciences location."
        )

    def shutdown(self):
        """shutdown deletes build packages if the user didn't request to clean
        up the cache. At this point we've already cancelled running jobs.
        """
        # Delete build source packages only if user regooglquested no cache
        if self._save_storage_cache:
            logger.debug("Requested to save workflow sources, skipping cleanup.")
        else:
            for package in self._build_packages:
                blob = self.bucket.blob(package)
                if blob.exists():
                    logger.debug("Deleting blob %s" % package)
                    blob.delete()

        # perform additional steps on shutdown if necessary
        super().shutdown()

    def cancel(self):
        """cancel execution, usually by way of control+c. Cleanup is done in
        shutdown (deleting cached workdirs in Google Cloud Storage
        """
        import googleapiclient

        # projects.locations.operations/cancel
        operations = self._api.projects().locations().operations()

        for job in self.active_jobs:
            request = operations.cancel(name=job.jobname)
            logger.debug("Cancelling operation {}".format(job.jobid))
            try:
                self._retry_request(request)
            except (Exception, BaseException, googleapiclient.errors.HttpError):
                continue

        self.shutdown()

    def get_available_machine_types(self):
        """Using the regions available at self.regions, use the GCP API
        to retrieve a lookup dictionary of all available machine types.
        """
        # Regular expression to determine if zone in region
        regexp = "^(%s)" % "|".join(self.regions)

        # Retrieve zones, filter down to selected regions
        zones = self._retry_request(
            self._compute_cli.zones().list(project=self.project)
        )
        zones = [z for z in zones["items"] if re.search(regexp, z["name"])]

        # Retrieve machine types available across zones
        # https://cloud.google.com/compute/docs/regions-zones/
        lookup = {}
        for zone in zones:
            request = self._compute_cli.machineTypes().list(
                project=self.project, zone=zone["name"]
            )
            lookup[zone["name"]] = self._retry_request(request)["items"]

        # Only keep those that are shared, use last zone as a base
        machine_types = {mt["name"]: mt for mt in lookup[zone["name"]]}
        del lookup[zone["name"]]

        # Update final list based on the remaining
        to_remove = set()
        for zone, types in lookup.items():
            names = [x["name"] for x in types]
            names = [name for name in names if "micro" not in name]
            names = [name for name in names if not re.search("^(e2|m1)", name)]
            for machine_type in list(machine_types.keys()):
                if machine_type not in names:
                    to_remove.add(machine_type)

        for machine_type in to_remove:
            del machine_types[machine_type]
        return machine_types

    def _add_gpu(self, gpu_count):
        """Add a number of NVIDIA gpus to the current executor. This works
        by way of adding nvidia_gpu to the job default resources, and also
        changing the default machine type prefix to be n1, which is
        the currently only supported instance type for using GPUs for LHS.
        """
        if not gpu_count or gpu_count == 0:
            return

        logger.debug(
            "found resource request for {} GPUs. This will limit to n1 "
            "instance types.".format(gpu_count)
        )
        self.workflow.default_resources.parsed["nvidia_gpu"] = gpu_count
        self.workflow.default_resources.args.append("nvidia_gpu=%s" % gpu_count)
        self._machine_type_prefix = self._machine_type_prefix or ""
        if not self._machine_type_prefix.startswith("n1"):
            self._machine_type_prefix = "n1"

    def _set_preemptible_rules(self, preemption_default=None, preemptible_rules=None):
        """define a lookup dictionary for preemptible instance retries, which
        is supported by the Google Life Science API. The user can set a default
        for all steps, specify per step, or define a default for all steps
        that aren't individually customized.
        """
        self.preemptible_rules = {}

        # If a default is defined, we apply it to all the rules
        if preemption_default is not None:
            self.preemptible_rules = {
                rule.name: preemption_default for rule in self.workflow.rules
            }

        # Now update custom set rules
        if preemptible_rules is not None:
            for rule in preemptible_rules:
                rule_name, restart_times = rule.strip().split("=")
                self.preemptible_rules[rule_name] = int(restart_times)

        # Ensure we set the number of restart times for each rule
        for rule_name, restart_times in self.preemptible_rules.items():
            rule = self.workflow.get_rule(rule_name)
            rule.restart_times = restart_times

    def _generate_job_resources(self, job):
        """given a particular job, generate the resources that it needs,
        including default regions and the virtual machine configuration
        """
        # Right now, do a best effort mapping of resources to instance types
        cores = job.resources.get("_cores", 1)
        mem_mb = job.resources.get("mem_mb", 15360)

        # IOPS performance proportional to disk size
        disk_mb = job.resources.get("disk_mb", 512000)

        # Convert mb to gb
        disk_gb = math.ceil(disk_mb / 1024)

        # Look for if the user wants an nvidia gpu
        gpu_count = job.resources.get("nvidia_gpu") or job.resources.get("gpu")
        gpu_model = job.resources.get("gpu_model")

        # If a gpu model is specified without a count, we assume 1
        if gpu_model and not gpu_count:
            gpu_count = 1

        # Update default resources using decided memory and disk
        self.workflow.default_resources = self.default_resources
        self.workflow.default_resources.args = [
            "mem_mb=%s" % mem_mb,
            "disk_mb=%s" % disk_mb,
        ]
        self.workflow.default_resources.parsed["mem_mb"] = mem_mb
        self.workflow.default_resources.parsed["disk_mb"] = disk_mb

        # Job resource specification can be overridden by gpu preferences
        self.machine_type_prefix = job.resources.get("machine_type")

        # If gpu wanted, limit to N1 general family, and update arguments
        if gpu_count:
            self._add_gpu(gpu_count)

        machine_types = self.get_available_machine_types()

        # Alert the user of machine_types available before filtering
        # https://cloud.google.com/compute/docs/machine-types
        logger.debug(
            "found {} machine types across regions {} before filtering "
            "to increase selection, define fewer regions".format(
                len(machine_types), self.regions
            )
        )

        # First pass - eliminate anything that too low in cpu/memory
        keepers = dict()

        # Also keep track of max cpus and memory, in case none available
        max_cpu = 1
        max_mem = 15360

        for name, machine_type in machine_types.items():
            max_cpu = max(max_cpu, machine_type["guestCpus"])
            max_mem = max(max_mem, machine_type["memoryMb"])
            if machine_type["guestCpus"] < cores or machine_type["memoryMb"] < mem_mb:
                continue
            keepers[name] = machine_type

        # If a prefix is set, filter down to it
        if self.machine_type_prefix:
            machine_types = keepers
            keepers = dict()
            for name, machine_type in machine_types.items():
                if name.startswith(self.machine_type_prefix):
                    keepers[name] = machine_type

        # If we don't have any contenders, workflow error
        if not keepers:
            if self.machine_type_prefix:
                raise WorkflowError(
                    "Machine prefix {prefix} is too strict, or the resources cannot "
                    " be satisfied, so there are no options "
                    "available.".format(prefix=self.machine_type_prefix)
                )
            else:
                raise WorkflowError(
                    "You requested {requestMemory} MB memory, {requestCpu} cores. "
                    "The maximum available are {availableMemory} MB memory and "
                    "{availableCpu} cores. These resources cannot be satisfied. "
                    "Please consider reducing the resource requirements of the "
                    "corresponding rule.".format(
                        requestMemory=mem_mb,
                        requestCpu=cores,
                        availableCpu=max_cpu,
                        availableMemory=max_mem,
                    )
                )

        # Now find (quasi) minimal to satisfy constraints
        machine_types = keepers

        # Select the first as the "smallest"
        smallest = list(machine_types.keys())[0]
        min_cores = machine_types[smallest]["guestCpus"]
        min_mem = machine_types[smallest]["memoryMb"]

        for name, machine_type in machine_types.items():
            if (
                machine_type["guestCpus"] < min_cores
                and machine_type["memoryMb"] < min_mem
            ):
                smallest = name
                min_cores = machine_type["guestCpus"]
                min_mem = machine_type["memoryMb"]

        selected = machine_types[smallest]
        logger.debug(
            "Selected machine type {}:{}".format(smallest, selected["description"])
        )

        virtual_machine = {
            "machineType": smallest,
            "labels": {"app": "snakemake"},
            "bootDiskSizeGb": disk_gb,
            "preemptible": job.rule.name in self.preemptible_rules,
        }

        # If the user wants gpus, add accelerators here
        if gpu_count:
            accelerator = self._get_accelerator(
                gpu_count, zone=selected["zone"], gpu_model=gpu_model
            )
            virtual_machine["accelerators"] = [
                {"type": accelerator["name"], "count": gpu_count}
            ]

        resources = {"regions": self.regions, "virtualMachine": virtual_machine}
        return resources

    def _get_accelerator(self, gpu_count, zone, gpu_model=None):
        """Get an appropriate accelerator for a GPU given a zone selection.
        Currently Google offers NVIDIA Tesla T4 (likely the best),
        NVIDIA P100, and the same T4 for a graphical workstation. Since
        this isn't a graphical workstation use case, we choose the
        accelerator that has >= to the maximumCardsPerInstace
        """
        if not gpu_count or gpu_count == 0:
            return

        accelerators = self._retry_request(
            self._compute_cli.acceleratorTypes().list(project=self.project, zone=zone)
        )

        # Filter down to those with greater than or equal to needed gpus
        keepers = {}
        for accelerator in accelerators.get("items", []):

            # Eliminate virtual workstations (vws) and models that don't match user preference
            if (gpu_model and accelerator["name"] != gpu_model) or accelerator[
                "name"
            ].endswith("vws"):
                continue

            if accelerator["maximumCardsPerInstance"] >= gpu_count:
                keepers[accelerator["name"]] = accelerator

        # If no matches available, exit early
        if not keepers:
            if gpu_model:
                raise WorkflowError(
                    "An accelerator in zone {zone} with model {model} cannot "
                    " be satisfied, so there are no options "
                    "available.".format(zone=zone, model=gpu_model)
                )
            else:
                raise WorkflowError(
                    "An accelerator in zone {zone} cannot be satisifed, so "
                    "there are no options available.".format(zone=zone)
                )

        # Find smallest (in future the user might have preference for the type)
        smallest = list(keepers.keys())[0]
        max_gpu = keepers[smallest]["maximumCardsPerInstance"]

        # This should usually return P-100, which would be preference (cheapest)
        for name, accelerator in keepers.items():
            if accelerator["maximumCardsPerInstance"] < max_gpu:
                smallest = name
                max_gpu = accelerator["maximumCardsPerInstance"]

        return keepers[smallest]

    def _set_snakefile(self):
        """The snakefile must be a relative path, which should be derived
        from the self.workflow.snakefile.
        """
        assert os.path.exists(self.workflow.snakefile)
        self.snakefile = self.workflow.snakefile.replace(self.workdir, "").strip(os.sep)

    def _set_workflow_sources(self):
        """We only add files from the working directory that are config related
        (e.g., the Snakefile or a config.yml equivalent), or checked into git.
        """
        self.workflow_sources = []

        for wfs in self.workflow.get_sources():
            if os.path.isdir(wfs):
                for (dirpath, dirnames, filenames) in os.walk(wfs):
                    self.workflow_sources.extend(
                        [
                            self.workflow.check_source_sizes(os.path.join(dirpath, f))
                            for f in filenames
                        ]
                    )
            else:
                self.workflow_sources.append(
                    self.workflow.check_source_sizes(os.path.abspath(wfs))
                )

    def _generate_build_source_package(self):
        """in order for the instance to access the working directory in storage,
        we need to upload it. This file is cleaned up at the end of the run.
        We do this, and then obtain from the instance and extract.
        """
        # Workflow sources for cloud executor must all be under same workdir root
        for filename in self.workflow_sources:
            if self.workdir not in filename:
                raise WorkflowError(
                    "All source files must be present in the working directory, "
                    "{workdir} to be uploaded to a build package that respects "
                    "relative paths, but {filename} was found outside of this "
                    "directory. Please set your working directory accordingly, "
                    "and the path of your Snakefile to be relative to it.".format(
                        workdir=self.workdir, filename=filename
                    )
                )

        # We will generate a tar.gz package, renamed by hash
        tmpname = next(tempfile._get_candidate_names())
        targz = os.path.join(tempfile.gettempdir(), "snakemake-%s.tar.gz" % tmpname)
        tar = tarfile.open(targz, "w:gz")

        # Add all workflow_sources files
        for filename in self.workflow_sources:
            arcname = filename.replace(self.workdir + os.path.sep, "")
            tar.add(filename, arcname=arcname)

        tar.close()

        # Rename based on hash, in case user wants to save cache
        sha256 = get_file_hash(targz)
        hash_tar = os.path.join(
            self.workflow.persistence.aux_path, "workdir-{}.tar.gz".format(sha256)
        )

        # Only copy if we don't have it yet, clean up if we do
        if not os.path.exists(hash_tar):
            shutil.move(targz, hash_tar)
        else:
            os.remove(targz)

        # We will clean these all up at shutdown
        self._build_packages.add(hash_tar)

        return hash_tar

    def _upload_build_source_package(self, targz):
        """given a .tar.gz created for a workflow, upload it to source/cache
        of Google storage, only if the blob doesn't already exist.
        """
        # Upload to temporary storage, only if doesn't exist
        self.pipeline_package = "source/cache/%s" % os.path.basename(targz)
        blob = self.bucket.blob(self.pipeline_package)
        logger.debug("build-package=%s" % self.pipeline_package)
        if not blob.exists():
            blob.upload_from_filename(targz, content_type="application/gzip")

    def _generate_log_action(self, job):
        """generate an action to save the pipeline logs to storage."""
        # script should be changed to this when added to version control!
        # https://raw.githubusercontent.com/snakemake/snakemake/master/snakemake/executors/google_lifesciences_helper.py

        # Save logs from /google/logs/output to source/logs in bucket
        commands = [
            "/bin/bash",
            "-c",
            "wget -O /gls.py https://raw.githubusercontent.com/snakemake/snakemake/master/snakemake/executors/google_lifesciences_helper.py && chmod +x /gls.py && source activate snakemake || true && python /gls.py save %s /google/logs %s/%s"
            % (self.bucket.name, self.gs_logs, job.name),
        ]

        # Always run the action to generate log output
        action = {
            "containerName": "snakelog-{}-{}".format(job.name, job.jobid),
            "imageUri": self.container_image,
            "commands": commands,
            "labels": self._generate_pipeline_labels(job),
            "alwaysRun": True,
        }
        return action

    def _generate_job_action(self, job):
        """generate a single action to execute the job."""
        # Derive the entrypoint command, the same content that might be written by self.get_jobscript(job)
        use_threads = "--force-use-threads" if not job.is_group() else ""

        exec_job = self.format_job(
            self.exec_job, job, _quote_all=True, use_threads=use_threads
        )

        # Now that we've parsed the job resource requirements, add to exec
        exec_job += self.get_default_resources_args()

        # script should be changed to this when added to version control!
        # https://raw.githubusercontent.com/snakemake/snakemake/master/snakemake/executors/google_lifesciences_helper.py
        # The full command to download the archive, extract, and run
        # For snakemake bases, we must activate the conda environment, but
        # for custom images we must allow this to fail (hence || true)
        commands = [
            "/bin/bash",
            "-c",
            "mkdir -p /workdir && cd /workdir && wget -O /download.py https://raw.githubusercontent.com/snakemake/snakemake/master/snakemake/executors/google_lifesciences_helper.py && chmod +x /download.py && source activate snakemake || true && python /download.py download %s %s /tmp/workdir.tar.gz && tar -xzvf /tmp/workdir.tar.gz && %s"
            % (self.bucket.name, self.pipeline_package, exec_job),
        ]

        # We are only generating one action, one job per run
        # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.pipelines/run#Action
        action = {
            "containerName": "snakejob-{}-{}".format(job.name, job.jobid),
            "imageUri": self.container_image,
            "commands": commands,
            "environment": self._generate_environment(),
            "labels": self._generate_pipeline_labels(job),
        }
        return action

    def _get_jobname(self, job):
        # Use a dummy job name (human readable and also namespaced)
        return "snakejob-%s-%s-%s" % (self.run_namespace, job.name, job.jobid)

    def _generate_pipeline_labels(self, job):
        """generate basic labels to identify the job, namespace, and that
        snakemake is running the show!
        """
        jobname = self._get_jobname(job)
        labels = {"name": jobname, "app": "snakemake"}
        return labels

    def _generate_environment(self):
        """loop through envvars (keys to host environment) and add
        any that are requested for the container environment.
        """
        envvars = {}
        for key in self.envvars:
            try:
                envvars[key] = os.environ[key]
            except KeyError:
                continue

        # Warn the user that we cannot support secrets
        if envvars:
            logger.warning("This API does not support environment secrets.")
        return envvars

    def _generate_pipeline(self, job):
        """based on the job details, generate a google Pipeline object
        to pass to pipelines.run. This includes actions, resources,
        environment, and timeout.
        """
        # Generate actions (one per job step) and log saving action (runs no matter what) and resources
        resources = self._generate_job_resources(job)
        action = self._generate_job_action(job)
        log_action = self._generate_log_action(job)

        pipeline = {
            # Ordered list of actions to execute
            "actions": [action, log_action],
            # resources required for execution
            "resources": resources,
            # Technical question - difference between resource and action environment
            # For now we will set them to be the same.
            "environment": self._generate_environment(),
        }

        # "timeout": string in seconds (3.5s) is not included (defaults to 7 days)
        return pipeline

    def run(self, job, callback=None, submit_callback=None, error_callback=None):

        super()._run(job)

        # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.pipelines
        pipelines = self._api.projects().locations().pipelines()

        # pipelines.run
        # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.pipelines/run

        labels = self._generate_pipeline_labels(job)
        pipeline = self._generate_pipeline(job)

        # The body of the request is a Pipeline and labels
        body = {"pipeline": pipeline, "labels": labels}

        # capabilities - this won't currently work (Singularity in Docker)
        # We either need to add CAPS or run in privileged mode (ehh)
        if job.needs_singularity and self.workflow.use_singularity:
            logger.warning(
                "Singularity requires additional capabilities that "
                "aren't yet supported for standard Docker runs."
            )

        # location looks like: "projects/<project>/locations/<location>"
        operation = pipelines.run(parent=self.location, body=body)

        # 403 will result if no permission to use pipelines or project
        result = self._retry_request(operation)

        # The jobid is the last number of the full name
        jobid = result["name"].split("/")[-1]

        # Give some logging for how to get status
        logger.info(
            "Get status with:\n"
            "gcloud config set project {project}\n"
            "gcloud beta lifesciences operations describe {location}/operations/{jobid}\n"
            "gcloud beta lifesciences operations list\n"
            "Logs will be saved to: {bucket}/{logdir}\n".format(
                project=self.project,
                jobid=jobid,
                location=self.location,
                bucket=self.bucket.name,
                logdir=self.gs_logs,
            )
        )

        self.active_jobs.append(
            GoogleLifeSciencesJob(job, result["name"], jobid, callback, error_callback)
        )

    def _job_was_successful(self, status):
        """based on a status response (a [pipeline].projects.locations.operations.get
        debug print the list of events, return True if all return codes 0
        and False otherwise (indication of failure). In that a nonzero exit
        status is found, we also debug print it for the user.
        """
        success = True

        # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/Event
        for event in status["metadata"]["events"]:

            logger.debug(event["description"])

            # Does it always result in fail for other failure reasons?
            if "failed" in event:
                success = False
                action = event.get("failed")
                logger.debug("{}: {}".format(action["code"], action["cause"]))

            elif "unexpectedExitStatus" in event:
                action = event.get("unexpectedExitStatus")

                if action["exitStatus"] != 0:
                    success = False

                    # Provide reason for the failure (desc includes exit code)
                    msg = "%s" % event["description"]
                    if "stderr" in action:
                        msg += ": %s" % action["stderr"]
                        logger.debug(msg)

        return success

    def _retry_request(self, request, timeout=2, attempts=3):
        """The Google Python API client frequently has BrokenPipe errors. This
        function takes a request, and executes it up to number of retry,
        each time with a 2* increase in timeout.

        Parameters
        ==========
        request: the Google Cloud request that needs to be executed
        timeout: time to sleep (in seconds) before trying again
        attempts: remaining attempts, throw error when hit 0
        """
        import googleapiclient

        try:
            return request.execute()
        except BrokenPipeError as ex:
            if attempts > 0:
                time.sleep(timeout)
                return self._retry_request(request, timeout * 2, attempts - 1)
            raise ex
        except googleapiclient.errors.HttpError as ex:
            log_verbose_traceback(ex)
            raise ex
        except Exception as ex:
            log_verbose_traceback(ex)
            raise ex

    def _wait_for_jobs(self):
        """wait for jobs to complete. This means requesting their status,
        and then marking them as finished when a "done" parameter
        shows up. Even for finished jobs, the status should still return
        """
        import googleapiclient

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

                # use self.status_rate_limiter to avoid too many API calls.
                with self.status_rate_limiter:

                    # https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.operations/get
                    # Get status from projects.locations.operations/get
                    operations = self._api.projects().locations().operations()
                    request = operations.get(name=j.jobname)
                    logger.debug("Checking status for operation {}".format(j.jobid))

                    try:
                        status = self._retry_request(request)
                    except googleapiclient.errors.HttpError as ex:

                        # Operation name not found, even finished should be found
                        if ex.status == 404:
                            j.error_callback(j.job)
                            continue

                        # Unpredictable server (500) error
                        elif ex.status == 500:
                            logger.error(ex["content"].decode("utf-8"))
                            j.error_callback(j.job)

                    except WorkflowError as ex:
                        print_exception(ex, self.workflow.linemaps)
                        j.error_callback(j.job)
                        continue

                    # The operation is done
                    if status.get("done", False) == True:

                        # Derive success/failure from status codes (prints too)
                        if self._job_was_successful(status):
                            j.callback(j.job)
                        else:
                            self.print_job_error(j.job, jobid=j.jobid)
                            j.error_callback(j.job)

                    # The operation is still running
                    else:
                        still_running.append(j)

            with self.lock:
                self.active_jobs.extend(still_running)
            sleep()
