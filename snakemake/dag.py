__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import asyncio
from builtins import ExceptionGroup
import datetime
import hashlib
import html
import os
import shutil
import subprocess
import tarfile
import textwrap
import time
from typing import Iterable, List, Optional, Set, Union
import uuid
import subprocess
from collections import Counter, defaultdict, deque, namedtuple
from functools import partial
from itertools import chain, filterfalse, groupby
from operator import attrgetter
from pathlib import Path
from snakemake.settings.types import DeploymentMethod

from snakemake_interface_executor_plugins.dag import DAGExecutorInterface
from snakemake_interface_report_plugins.interfaces import DAGReportInterface
from snakemake_interface_storage_plugins.storage_object import StorageObjectTouch

from snakemake import workflow
from snakemake import workflow as _workflow
from snakemake.common import (
    ON_WINDOWS,
    group_into_chunks,
    is_local_file,
)
from snakemake.settings.types import RerunTrigger
from snakemake.deployment import singularity
from snakemake.exceptions import (
    AmbiguousRuleException,
    ChildIOException,
    CyclicGraphException,
    ImproperOutputException,
    IncompleteFilesException,
    InputFunctionException,
    MissingInputException,
    MissingOutputException,
    MissingRuleException,
    PeriodicWildcardError,
    RemoteFileException,
    WildcardError,
    WorkflowError,
)
from snakemake.io import (
    _IOFile,
    PeriodicityDetector,
    get_flag_value,
    is_callable,
    is_flagged,
    wait_for_files,
)
from snakemake.jobs import (
    AbstractJob,
    GroupJob,
    GroupJobFactory,
    Job,
    JobFactory,
    Reason,
)
from snakemake.settings.types import SharedFSUsage
from snakemake.logging import logger
from snakemake.output_index import OutputIndex
from snakemake.sourcecache import LocalSourceFile, SourceFile
from snakemake.settings.types import ChangeType, Batch

PotentialDependency = namedtuple("PotentialDependency", ["file", "jobs", "known"])


def toposort(graph):
    from graphlib import TopologicalSorter

    sorter = TopologicalSorter(graph)
    sorter.prepare()
    sorted = list()
    while sorter.is_active():
        ready = set()
        for task in sorter.get_ready():
            ready.add(task)
            sorter.done(task)
        sorted.append(ready)
    return sorted


class DAG(DAGExecutorInterface, DAGReportInterface):
    """Directed acyclic graph of jobs."""

    def __init__(
        self,
        workflow,
        rules=None,
        targetfiles: Set[str] = None,
        targetrules=None,
        forceall=False,
        forcerules=None,
        forcefiles=None,
        priorityfiles=None,
        priorityrules=None,
        untilfiles=None,
        untilrules=None,
        omitfiles=None,
        omitrules=None,
        ignore_incomplete=False,
    ):
        self._queue_input_jobs = None
        self._dependencies = defaultdict(partial(defaultdict, set))
        self.depending = defaultdict(partial(defaultdict, set))
        self._needrun = set()
        self._checkpoint_jobs = set()
        self._priority = dict()
        self._reason = defaultdict(Reason)
        self._finished = set()
        self._has_unfinished_queue_input_jobs = None
        self._len = 0
        self.workflow: _workflow.Workflow = workflow
        self.rules = set(rules)
        self.targetfiles = targetfiles
        self.targetrules = targetrules
        self.target_jobs_rules = {
            spec.rulename for spec in self.workflow.dag_settings.target_jobs
        }
        self.priorityfiles = priorityfiles
        self.priorityrules = priorityrules
        self.targetjobs = set()
        self.prioritytargetjobs = set()
        self._ready_jobs = set()
        self._jobid = dict()
        self.job_cache = dict()
        self.conda_envs = dict()
        self.container_imgs = dict()
        self._progress = 0
        self._group = dict()
        self._n_until_ready = defaultdict(int)
        self._running = set()
        self._jobs_with_finished_queue_input = set()
        self._storage_input_jobs = defaultdict(list)

        self.job_factory = JobFactory()
        self.group_job_factory = GroupJobFactory()

        self.forcerules = set()
        self.forcefiles = set()
        self.untilrules = set()
        self.untilfiles = set()
        self.omitrules = set()
        self.omitfiles = set()
        if forceall:
            self.forcerules.update(self.rules)
        elif forcerules:
            self.forcerules.update(forcerules)
        if forcefiles:
            self.forcefiles.update(forcefiles)
        if untilrules:
            self.untilrules.update(set(rule.name for rule in untilrules))
        if untilfiles:
            self.untilfiles.update(untilfiles)
        if omitrules:
            self.omitrules.update(set(rule.name for rule in omitrules))
        if omitfiles:
            self.omitfiles.update(omitfiles)

        self.omitforce = set()

        if self.batch is not None and not self.batch.is_final:
            # Since not all input files of a batching rule are considered, we cannot run
            # beyond that rule.
            # For the final batch, we do not need to omit anything.
            self.omitrules.add(self.batch.rulename)

        self.ignore_incomplete = ignore_incomplete

        self.periodic_wildcard_detector = PeriodicityDetector()

        self.update_output_index()

    @property
    def dependencies(self):
        return self._dependencies

    @property
    def batch(self):
        return self.workflow.dag_settings.batch

    async def init(self, progress=False):
        """Initialise the DAG."""
        for job in [await self.rule2job(rule) for rule in self.targetrules]:
            job = await self.update([job], progress=progress, create_inventory=True)
            self.targetjobs.add(job)

        for file in self.targetfiles:
            job = await self.update(
                await self.file2jobs(file),
                file=file,
                progress=progress,
                create_inventory=True,
            )
            self.targetjobs.add(job)

        for spec in self.workflow.dag_settings.target_jobs:
            job = await self.update(
                [
                    await self.new_job(
                        self.workflow.get_rule(spec.rulename),
                        wildcards_dict=spec.wildcards_dict,
                    )
                ],
                progress=progress,
                create_inventory=True,
            )
            self.targetjobs.add(job)
            self.forcefiles.update(job.output)

        self.cleanup()

        await self.check_incomplete()

        self.update_container_imgs()
        self.update_conda_envs()

        await self.update_needrun(create_inventory=True)
        self.set_until_jobs()
        self.delete_omitfrom_jobs()
        self.update_jobids()

        self.check_directory_outputs()

        # check if remaining jobs are valid
        for i, job in enumerate(self.jobs):
            job.is_valid()

    def get_unneeded_temp_files(self, job: AbstractJob) -> Iterable[str]:
        if isinstance(job, GroupJob):
            for j in job:
                yield from self.get_unneeded_temp_files(j)
        else:
            for f in job.output:
                if is_flagged(f, "temp") and not self.is_needed_tempfile(job, f):
                    yield f

    def check_directory_outputs(self):
        """Check that no output file is contained in a directory output of the same or another rule."""
        outputs = sorted(
            {(os.path.abspath(f), job) for job in self.jobs for f in job.output}
        )
        for i in range(len(outputs) - 1):
            (a, job_a), (b, job_b) = outputs[i : i + 2]
            try:
                common = os.path.commonpath([a, b])
            except ValueError:
                # commonpath raises error if windows drives are different.
                continue
            if a != b and common == os.path.commonpath([a]) and job_a != job_b:
                raise ChildIOException(parent=outputs[i], child=outputs[i + 1])

    @property
    def checkpoint_jobs(self):
        return self._checkpoint_jobs

    @property
    def finished_checkpoint_jobs(self):
        for job in self.finished_jobs:
            if job.is_checkpoint:
                yield job

    def update_checkpoint_outputs(self):
        workflow.checkpoints.future_output = set(
            f for job in self.checkpoint_jobs for f in job.output
        )
        workflow.checkpoints.created_output = set(
            f for job in self.finished_checkpoint_jobs for f in job.output
        )

    def update_jobids(self):
        for job in self.jobs:
            if job not in self._jobid:
                self._jobid[job] = len(self._jobid)

    def cleanup_workdir(self):
        for job in self.jobs:
            if not self.is_edit_notebook_job(job):
                for io_dir in set(
                    os.path.dirname(io_file)
                    for io_file in chain(job.output, job.input)
                    if not os.path.exists(io_file)
                ):
                    if os.path.exists(io_dir):
                        # check for empty dir
                        with os.scandir(io_dir) as i:
                            if next(i, None) is None:
                                os.removedirs(io_dir)

    def cleanup(self):
        self.job_cache.clear()
        final_jobs = set(self.bfs(self._dependencies, *self.targetjobs))
        todelete = [job for job in self._dependencies if job not in final_jobs]
        for job in todelete:
            try:
                self._needrun.remove(job)
            except KeyError:
                pass

            # delete all pointers from dependencies to this job
            for dep in self._dependencies[job]:
                try:
                    del self.depending[dep][job]
                except KeyError:
                    # In case the pointer has been deleted before or
                    # never created, we can simply continue.
                    pass

            # delete all dependencies
            del self._dependencies[job]
            try:
                # delete all pointers to downstream dependencies
                del self.depending[job]
            except KeyError:
                pass

    def update_conda_envs(self):
        # First deduplicate based on job.conda_env_spec
        env_set = {
            (job.conda_env_spec, job.container_img_url)
            for job in self.jobs
            if job.conda_env_spec
            and (
                job.is_local
                or SharedFSUsage.SOFTWARE_DEPLOYMENT
                in self.workflow.storage_settings.shared_fs_usage
                or (
                    self.workflow.remote_exec
                    and SharedFSUsage.SOFTWARE_DEPLOYMENT
                    not in self.workflow.storage_settings.shared_fs_usage
                )
            )
        }

        # Then based on md5sum values
        for env_spec, simg_url in env_set:
            simg = None
            if simg_url and (
                DeploymentMethod.APPTAINER
                in self.workflow.deployment_settings.deployment_method
            ):
                assert (
                    simg_url in self.container_imgs
                ), "bug: must first pull singularity images"
                simg = self.container_imgs[simg_url]
            key = (env_spec, simg_url)
            if key not in self.conda_envs:
                env = env_spec.get_conda_env(
                    self.workflow,
                    container_img=simg,
                    cleanup=self.workflow.deployment_settings.conda_cleanup_pkgs,
                )
                self.conda_envs[key] = env

    async def retrieve_storage_inputs(self, jobs=None, also_missing_internal=False):
        shared_local_copies = (
            SharedFSUsage.STORAGE_LOCAL_COPIES
            in self.workflow.storage_settings.shared_fs_usage
        )
        if jobs is None:
            if (self.workflow.is_main_process and shared_local_copies) or (
                self.workflow.remote_exec and not shared_local_copies
            ):
                jobs = self.needrun_jobs()
            else:
                jobs = []

        to_retrieve = {
            f
            for job in jobs
            for f in job.input
            if f.is_storage
            and (
                (also_missing_internal and not shared_local_copies)
                or self.is_external_input(f, job, not_needrun_is_external=True)
            )
            and not job.is_norun
        }

        if to_retrieve:
            try:
                async with asyncio.TaskGroup() as tg:
                    for f in to_retrieve:
                        logger.info(f"Retrieving {f} from storage.")
                        tg.create_task(f.retrieve_from_storage())
            except ExceptionGroup as e:
                raise WorkflowError("Failed to retrieve input from storage.", e)

    def update_storage_inputs(self):
        self._storage_input_jobs.clear()
        for job in self.needrun_jobs():
            for f in job.input:
                if f.is_storage:
                    self._storage_input_jobs[f].append(job)

    async def store_storage_outputs(self):
        if self.workflow.remote_exec:
            logger.info("Storing output in storage.")
            try:
                async with asyncio.TaskGroup() as tg:
                    for job in self.needrun_jobs(exclude_finished=False):
                        benchmark = [job.benchmark] if job.benchmark else []

                        async def tostore(f):
                            return (
                                f.is_storage
                                and f
                                not in self.workflow.storage_settings.unneeded_temp_files
                                and await f.exists_local()
                            )

                        if self.finished(job):
                            for f in chain(job.output, benchmark):
                                if await tostore(f):
                                    tg.create_task(f.store_in_storage())
                        for f in job.log:
                            if await tostore(f):
                                tg.create_task(f.store_in_storage())
            except ExceptionGroup as e:
                raise WorkflowError("Failed to store output in storage.", e)

    async def cleanup_storage_objects(self):
        shared_local_copies = (
            SharedFSUsage.STORAGE_LOCAL_COPIES
            in self.workflow.storage_settings.shared_fs_usage
        )
        cleaned = set()
        for job in self.jobs:
            if (
                self.workflow.is_main_process and (job.is_local or shared_local_copies)
            ) or (self.workflow.remote_exec and not shared_local_copies):
                async with asyncio.TaskGroup() as tg:
                    for f in chain(job.input, job.output):
                        if f.is_storage and f not in cleaned:
                            f.storage_object.cleanup()
                            tg.create_task(
                                f.remove(remove_non_empty_dir=True, only_local=True)
                            )
                            cleaned.add(f)

    async def sanitize_local_storage_copies(self):
        """Remove local copies of storage files that will be recreated in this run."""
        async with asyncio.TaskGroup() as tg:
            for job in self.needrun_jobs():
                if not self.finished(job):
                    for f in job.output:
                        if f.is_storage and await f.exists_local():
                            tg.create_task(
                                f.remove(remove_non_empty_dir=True, only_local=True)
                            )

    def create_conda_envs(self, dryrun=False, quiet=False):
        dryrun |= self.workflow.dryrun
        for env in self.conda_envs.values():
            if (not dryrun or not quiet) and not env.is_externally_managed:
                env.create(self.workflow.dryrun)

    def update_container_imgs(self):
        # First deduplicate based on job.conda_env_spec
        img_set = {
            (job.container_img_url, job.is_containerized)
            for job in self.jobs
            if job.container_img_url
        }

        for img_url, is_containerized in img_set:
            if img_url not in self.container_imgs:
                img = singularity.Image(img_url, self, is_containerized)
                self.container_imgs[img_url] = img

    def pull_container_imgs(self, quiet=False):
        for img in self.container_imgs.values():
            if not self.workflow.dryrun or not quiet:
                img.pull(self.workflow.dryrun)

    def update_output_index(self):
        """Update the OutputIndex."""
        self.output_index = OutputIndex(self.rules)

    async def check_incomplete(self):
        """Check if any output files are incomplete. This is done by looking up
        markers in the persistence module."""
        if not self.ignore_incomplete:
            incomplete = await self.incomplete_files()
            if incomplete:
                if self.workflow.dag_settings.force_incomplete:
                    self.forcefiles.update(incomplete)
                else:
                    raise IncompleteFilesException(incomplete)

    def incomplete_external_jobid(self, job) -> Optional[str]:
        """Return the external jobid of the job if it is marked as incomplete.

        Returns None, if job is not incomplete, or if no external jobid has been
        registered or if force_incomplete is True.
        """
        if self.workflow.dag_settings.force_incomplete:
            return None
        jobids = self.workflow.persistence.external_jobids(job)
        if len(jobids) == 1:
            return jobids[0]
        elif len(jobids) > 1:
            raise WorkflowError(
                "Multiple different external jobids registered "
                "for output files of incomplete job {} ({}). This job "
                "cannot be resumed. Execute Snakemake with --rerun-incomplete "
                "to fix this issue.".format(job.jobid, jobids)
            )

    def is_edit_notebook_job(self, job):
        return (
            self.workflow.execution_settings.edit_notebook
            and job.targetfile in self.targetfiles
        )

    def is_draft_notebook_job(self, job):
        return (
            self.workflow.execution_settings.edit_notebook
            and self.workflow.execution_settings.edit_notebook.draft_only
            and job.targetfile in self.targetfiles
        )

    def get_job_group(self, job):
        return self._group.get(job)

    @property
    def jobs(self):
        """All jobs in the DAG."""
        return self._dependencies.keys()

    def needrun_jobs(self, exclude_finished=True):
        """Jobs that need to be executed."""
        if exclude_finished:
            return filterfalse(self.finished, self._needrun)
        else:
            return iter(self._needrun)

    @property
    def local_needrun_jobs(self):
        """Iterate over all jobs that need to be run and are marked as local."""
        return filter(lambda job: job.is_local, self.needrun_jobs())

    @property
    def finished_jobs(self):
        """Iterate over all jobs that have been finished."""
        return filter(self.finished, self.jobs)

    @property
    def ready_jobs(self):
        """Jobs that are ready to execute."""
        return self._ready_jobs

    def needrun(self, job):
        """Return whether a given job needs to be executed."""
        return job in self._needrun

    def priority(self, job):
        """Return priority of given job."""
        return self._priority[job]

    def noneedrun_finished(self, job):
        """
        Return whether a given job is finished or was not
        required to run at all.
        """
        return not self.needrun(job) or self.finished(job)

    def reason(self, job):
        """Return the reason of the job execution."""
        return self._reason[job]

    def finished(self, job):
        """Return whether a job is finished."""
        return job in self._finished

    def requested_files(self, job):
        """Return the files a job requests."""
        return set(*self.depending[job].values())

    async def incomplete_files(self):
        """Yield incomplete files."""
        incomplete = list()
        for job in filterfalse(self.needrun, self.jobs):
            is_incomplete = await self.workflow.persistence.incomplete(job)
            if is_incomplete:
                for f in job.output:
                    incomplete.append(f)
        return incomplete

    @property
    def newversion_files(self):
        """Return list of files where the current version is newer than the
        recorded version.
        """
        return list(
            chain(
                *(
                    job.output
                    for job in filter(self.workflow.persistence.newversion, self.jobs)
                )
            )
        )

    async def handle_ensure(self, job, expanded_output):
        ensured_output = {
            f: get_flag_value(f, "ensure")
            for f in expanded_output
            if is_flagged(f, "ensure")
        }
        # handle non_empty
        empty_output = [
            f
            for f, ensure in ensured_output.items()
            if ensure["non_empty"] and (await f.size()) == 0
        ]
        if empty_output:
            raise WorkflowError(
                "Detected unexpected empty output files. "
                "Something went wrong in the rule without "
                "an error being reported:\n{}".format("\n".join(empty_output)),
                rule=job.rule,
            )

        # handle checksum
        async def is_not_same_checksum(f, checksum):
            if checksum is None:
                return False
            if is_callable(checksum):
                try:
                    checksum = checksum(job.wildcards)
                except Exception as e:
                    raise WorkflowError(
                        "Error calling checksum function provided to ensure marker.",
                        e,
                        rule=job.rule,
                    )
            return not await f.is_same_checksum(checksum, force=True)

        checksum_failed_output = [
            f
            for f, ensure in ensured_output.items()
            if await is_not_same_checksum(f, ensure.get("sha256"))
        ]
        if checksum_failed_output:
            raise WorkflowError(
                "Output files have checksums that differ from the expected ones "
                "defined in the workflow:\n{}".format(
                    "\n".join(checksum_failed_output)
                ),
                rule=job.rule,
            )

    def check_touch_compatible(self):
        def is_touchable(f):
            return not f.is_storage or isinstance(f.storage_object, StorageObjectTouch)

        if not all(is_touchable(f) for job in self.jobs for f in job.output):
            raise WorkflowError(
                "Touching output files is impossible. The workflow uses remote storage "
                "but the storage plugin does not support the touch operation. "
                "It might be possible to improve the storage plugin to support this "
                "operation. Consider checking the source code and contributing to the "
                "plugin."
            )

    async def check_and_touch_output(
        self,
        job,
        wait=3,
        ignore_missing_output=False,
        no_touch=False,
        wait_for_local=True,
    ):
        """Raise exception if output files of job are missing."""
        # do not touch output in storage. This is done before by the touch executor.
        expanded_output = [job.shadowed_path(path) for path in job.output]
        if job.benchmark:
            expanded_output.append(job.benchmark)

        if not ignore_missing_output:
            try:
                await wait_for_files(
                    expanded_output,
                    latency_wait=wait,
                    wait_for_local=wait_for_local,
                    ignore_pipe_or_service=True,
                )
            except IOError as e:
                raise MissingOutputException(
                    str(e), rule=job.rule, jobid=self.jobid(job)
                )

        def correctly_flagged_with_dir(f):
            """Check that files flagged as directories are in fact directories

            In ambiguous cases, such as when f is managed by a storage backend, or f
            doesn't exist and
            ignore_missing_output is true, always return True
            """
            if f.storage_object:
                return True
            if ignore_missing_output and not os.path.exists(f):
                return True
            return not (f.is_directory ^ os.path.isdir(f))

        for output_path in expanded_output:
            if not correctly_flagged_with_dir(output_path):
                raise ImproperOutputException(job, [output_path])

        # Handle ensure flags
        await self.handle_ensure(job, expanded_output)

        # It is possible, due to archive expansion or cluster clock skew, that
        # the files appear older than the input.  But we know they must be new,
        # so touch them to update timestamps. This also serves to touch outputs
        # when using the --touch flag.
        if not no_touch:
            for output_path in expanded_output:
                # This won't create normal files if missing, but will create
                # the flag file for directories.
                if await output_path.exists_local():
                    output_path.touch()

        if wait_for_local:
            await self.check_output_mtime(job, expanded_output)

    async def check_output_mtime(
        self, job: Job, expanded_output: List[_IOFile]
    ) -> None:
        # Check whether all output files are newer than the newest input file.
        # This catches problems with non-sychronous clocks in distributed execution.
        existing_input = [
            ((await f.mtime_uncached(skip_storage=True)).local(), f)
            for f in job.input
            if await f.exists_local()
        ]
        if existing_input:
            newest_input_mtime, newest_input_path = max(existing_input)
            for output_path in expanded_output:
                if not await output_path.exists_local():
                    # can happen for temp files in combination with touch executor
                    continue
                output_mtime = (
                    await output_path.mtime_uncached(skip_storage=True)
                ).local()
                if output_mtime < newest_input_mtime:
                    raise WorkflowError(
                        f"Output {output_path} has older modification time "
                        f"({datetime.datetime.fromtimestamp(output_mtime)}) "
                        f"than input {newest_input_path} "
                        f"({datetime.datetime.fromtimestamp(newest_input_mtime)}). "
                        "This could indicate a clock skew problem in your network "
                        "and would trigger a rerun of this job in the next "
                        "execution and should therefore be fixed on system level. "
                        f"System time: {datetime.datetime.now()}",
                        rule=job.rule,
                    )

    def unshadow_output(self, job, only_log=False):
        """Move files from shadow directory to real output paths."""
        if not job.shadow_dir or not job.output:
            return

        files = job.log if only_log else chain(job.output, job.log)

        for real_output in files:
            shadow_output = job.shadowed_path(real_output).file
            # Remake absolute symlinks as relative
            if os.path.islink(shadow_output):
                dest = os.readlink(shadow_output)
                if os.path.isabs(dest):
                    rel_dest = os.path.relpath(dest, job.shadow_dir)
                    os.remove(shadow_output)
                    os.symlink(rel_dest, shadow_output)

            if os.path.realpath(shadow_output) == os.path.realpath(real_output):
                continue
            shutil.move(shadow_output, real_output)
        shutil.rmtree(job.shadow_dir)

    def check_periodic_wildcards(self, job):
        """Raise an exception if a wildcard of the given job appears to be periodic,
        indicating a cyclic dependency."""
        for wildcard, value in job.wildcards_dict.items():
            periodic_substring = self.periodic_wildcard_detector.is_periodic(value)
            if periodic_substring is not None:
                raise PeriodicWildcardError(
                    "The value {} in wildcard {} is periodically repeated ({}). "
                    "This would lead to an infinite recursion. "
                    "To avoid this, e.g. restrict the wildcards in this rule to certain values.".format(
                        periodic_substring, wildcard, value
                    ),
                    rule=job.rule,
                )

    def handle_protected(self, job):
        """Write-protect output files that are marked with protected()."""
        for f in job.output:
            if f in job.protected_output:
                logger.info(f"Write-protecting output file {f}.")
                f.protect()

    def handle_touch(self, job):
        """Touches those output files that are marked for touching."""
        for f in job.output:
            if f in job.touch_output:
                f = job.shadowed_path(f)
                logger.info(f"Touching output file {f}.")
                f.touch_or_create()
                assert os.path.exists(f)

    def temp_input(self, job: Union[Job, GroupJob]):
        if job.is_group():
            skip = {f for j in job for f in j.output}
            jobs = job
        else:
            skip = set()
            jobs = [job]
        for job in jobs:
            for job_, files in self._dependencies[job].items():
                for f in files:
                    if f in job_.temp_output and f not in skip:
                        yield f

    async def temp_size(self, job):
        """Return the total size of temporary input files of the job.
        If none, return 0.
        """
        return sum([await f.size() for f in self.temp_input(job)])

    def is_needed_tempfile(self, job, tempfile):
        return any(
            tempfile in files
            for j, files in self.depending[job].items()
            if not self.finished(j) and self.needrun(j) and j != job
        )

    async def handle_temp(self, job):
        """Remove temp files if they are no longer needed. Update temp_mtimes."""
        if self.workflow.storage_settings.notemp:
            return

        if job.is_group():
            for j in job:
                await self.handle_temp(j)
            return

        is_temp = lambda f: is_flagged(f, "temp")

        def unneeded_files():
            # temp input
            for job_, files in self._dependencies[job].items():
                tempfiles = set(f for f in job_.output if is_temp(f))
                yield from filterfalse(
                    partial(self.is_needed_tempfile, job_), tempfiles & files
                )

            # temp output
            if not job.is_checkpoint and (
                job not in self.targetjobs
                or job.rule.name == self.workflow.default_target
            ):
                tempfiles = (
                    f for f in job.output if is_temp(f) and f not in self.targetfiles
                )
                yield from filterfalse(partial(self.is_needed_tempfile, job), tempfiles)

        for f in unneeded_files():
            if self.workflow.dryrun:
                logger.info(f"Would remove temporary output {f}")
            else:
                logger.info(f"Removing temporary output {f}.")
                await f.remove(remove_non_empty_dir=True)

    async def handle_log(self, job):
        for f in job.log:
            f = job.shadowed_path(f)
            if not await f.exists_local():
                # If log file was not created during job, create an empty one.
                f.touch_or_create()

    async def handle_storage(self, job, store_in_storage=True, store_only_log=False):
        """Remove local files if they are no longer needed and upload."""
        if store_in_storage and (
            self.workflow.remote_exec or self.workflow.is_main_process
        ):
            # handle output files
            files = []
            if not store_only_log:
                files = job.output
                if job.benchmark:
                    files = chain(files, (job.benchmark,))
            if job.log:
                files = chain(files, job.log)
            for f in files:
                if f.is_storage and not f.should_not_be_retrieved_from_storage:
                    await f.store_in_storage()
                    storage_mtime = (await f.mtime()).storage()
                    # immediately force local mtime to match storage,
                    # since conversions from S3 headers are not 100% reliable
                    # without this, newness comparisons may fail down the line
                    f.touch(times=(storage_mtime, storage_mtime))

                    if not await f.exists_in_storage():
                        raise RemoteFileException(
                            "The file upload was attempted, but it does not "
                            "exist in storage. Check that your credentials have "
                            "read AND write permissions."
                        )

        if not self.workflow.storage_settings.keep_storage_local:
            if not any(f.is_storage for f in job.input):
                return

            # handle input files
            needed = lambda job_, f: any(
                f in files
                for j, files in self.depending[job_].items()
                if not self.finished(j) and self.needrun(j) and j != job
            )

            async def unneeded_files():
                async def putative(f):
                    return (
                        f.is_storage
                        and not await f.protected()
                        and not f.should_keep_local
                    )

                generated_input = set()
                for job_, files in self._dependencies[job].items():
                    generated_input |= files
                    for f in files:
                        if await putative(f) and not needed(job_, f):
                            yield f
                for f, f_ in zip(job.output, job.rule.output):
                    if (
                        await putative(f)
                        and not needed(job, f)
                        and not f in self.targetfiles
                    ):
                        yield f
                for f in job.input:
                    if (
                        await putative(f)
                        and f not in generated_input
                        # all other jobs that depend on this file are finished
                        and all(
                            self.finished(job_)
                            for job_ in self._storage_input_jobs[f]
                            if job_ != job
                        )
                    ):
                        yield f

            async for f in unneeded_files():
                if await f.exists_local():
                    logger.info(f"Removing local copy of storage file: {f}")
                    await f.remove()

    def jobid(self, job):
        """Return job id of given job."""
        if job.is_group():
            return job.jobid
        else:
            return self._jobid[job]

    async def update(
        self,
        jobs,
        file=None,
        visited=None,
        known_producers=None,
        progress=False,
        create_inventory=False,
    ):
        """Update the DAG by adding given jobs and their dependencies."""
        if visited is None:
            visited = set()
        if known_producers is None:
            known_producers = dict()
        producers = []
        exceptions = list()
        cycles = list()

        # check if all potential producers are strictly ordered
        jobs = sorted(jobs, reverse=True)
        discarded_jobs = set()

        def is_strictly_higher_ordered(pivot_job):
            return all(
                pivot_job > job
                for job in jobs
                if job is not pivot_job and job not in discarded_jobs
            )

        for i, job in enumerate(jobs):
            logger.dag_debug(dict(status="candidate", job=job))
            if file in job.input:
                cycles.append(job)
                continue
            if job in visited:
                cycles.append(job)
                continue
            try:
                self.check_periodic_wildcards(job)
                await self.update_(
                    job,
                    visited=set(visited),
                    known_producers=known_producers,
                    progress=progress,
                    create_inventory=create_inventory,
                )
                producers.append(job)
                if is_strictly_higher_ordered(job):
                    # All other jobs are either discarded or strictly less given
                    # any defined ruleorder.
                    break
            except (
                MissingInputException,
                CyclicGraphException,
                PeriodicWildcardError,
                WorkflowError,
            ) as ex:
                exceptions.append(ex)
                discarded_jobs.add(job)
            except RecursionError as e:
                raise WorkflowError(
                    e,
                    (
                        "If building the DAG exceeds the recursion limit, "
                        "this is likely due to a cyclic dependency."
                        "E.g. you might have a sequence of rules that "
                        "can generate their own input. Try to make "
                        "the output files more specific. "
                        "A common pattern is to have different prefixes "
                        "in the output files of different rules."
                        + f"\nProblematic file pattern: {file}"
                        if file
                        else ""
                    ),
                )
        if not producers:
            if cycles:
                job = cycles[0]
                raise CyclicGraphException(job.rule, file, rule=job.rule)
            if len(exceptions) > 1:
                raise WorkflowError(*exceptions)
            elif len(exceptions) == 1:
                raise exceptions[0]

        n = len(self._dependencies)
        if progress and n % 1000 == 0 and n and self._progress != n:
            logger.info(f"Processed {n} potential jobs.")
            self._progress = n

        producers.sort(reverse=True)
        producer = producers[0]
        ambiguities = list(
            filter(lambda x: not x < producer and not producer < x, producers[1:])
        )
        if ambiguities and not self.workflow.execution_settings.ignore_ambiguity:
            raise AmbiguousRuleException(file, producer, ambiguities[0])
        logger.dag_debug(dict(status="selected", job=producer))
        if exceptions:
            logger.dag_debug(
                dict(
                    file=file,
                    msg="Producer found, hence exceptions are ignored.",
                    exception=WorkflowError(*exceptions),
                )
            )
        return producer

    async def update_(
        self,
        job,
        visited=None,
        known_producers=None,
        progress=False,
        create_inventory=False,
    ):
        """Update the DAG by adding the given job and its dependencies."""
        if job in self._dependencies:
            return
        if visited is None:
            visited = set()
        if known_producers is None:
            known_producers = dict()
        visited.add(job)
        dependencies = self._dependencies[job]
        potential_dependencies = [
            res
            async for res in self.collect_potential_dependencies(
                job, known_producers=known_producers
            )
        ]

        missing_input = set()
        producer = dict()
        exceptions = dict()
        for res in potential_dependencies:
            if create_inventory:
                # If possible, obtain inventory information starting from
                # given file and store it in the IOCache.
                # This should provide faster access to existence and mtime information
                # than querying file by file. If the file type does not support inventory
                # information, this call is a no-op.
                await res.file.inventory()

            if not res.jobs:
                # no producing job found
                if self.workflow.is_main_process and not await res.file.exists():
                    # file not found, hence missing input
                    missing_input.add(res.file)
                if not is_flagged(res.file, "before_update"):
                    # record info that there is no known producer, but only
                    # do that for files that are not flagged as before_update
                    # otherwise, we would risk a conflict with corresponding files
                    # that are flagged as 'update'.
                    known_producers[res.file] = None
                # file found, no problem
                continue

            if res.known:
                producer[res.file] = res.jobs[0]
            else:
                try:
                    selected_job = await self.update(
                        res.jobs,
                        file=res.file,
                        visited=visited,
                        known_producers=known_producers,
                        progress=progress,
                    )
                    producer[res.file] = selected_job
                except (
                    MissingInputException,
                    CyclicGraphException,
                    PeriodicWildcardError,
                    WorkflowError,
                ) as ex:
                    file_exists = await res.file.exists()
                    if not file_exists:
                        self.delete_job(job, recursive=False)  # delete job from tree
                        raise ex
                    else:
                        logger.dag_debug(
                            dict(
                                file=res.file,
                                msg="No producers found, but file is present on disk.",
                                exception=ex,
                            )
                        )
                        known_producers[res.file] = None

        for file, job_ in producer.items():
            dependencies[job_].add(file)
            self.depending[job_][job].add(file)

        if self.is_batch_rule(job.rule) and self.batch.is_final:
            # For the final batch, ensure that all input files from
            # previous batches are present on disk.
            if any([(f not in producer and not await f.exists()) for f in job.input]):
                raise WorkflowError(
                    "Unable to execute batch {} because not all previous batches "
                    "have been completed before or files have been deleted.".format(
                        self.batch
                    )
                )

        if missing_input:
            self.delete_job(job, recursive=False)  # delete job from tree
            raise MissingInputException(job, missing_input)

    async def update_needrun(self, create_inventory=False):
        """Update the information whether a job needs to be executed."""

        if create_inventory and self.workflow.is_main_process:
            # Concurrently collect mtimes of all existing files.
            await self.workflow.iocache.mtime_inventory(self.jobs)

        output_mintime = dict()

        async def update_output_mintime(job):
            try:
                return output_mintime[job]
            except KeyError:
                for job_ in chain([job], self.depending[job]):
                    try:
                        t = output_mintime[job_]
                    except KeyError:
                        t = await job_.output_mintime
                    if t is not None:
                        output_mintime[job] = t
                        return
                output_mintime[job] = None

        is_same_checksum_cache = dict()

        async def is_same_checksum(f, job):
            try:
                return is_same_checksum_cache[(f, job)]
            except KeyError:
                if not await f.is_checksum_eligible():
                    # no chance to compute checksum, cannot be assumed the same
                    is_same = False
                else:
                    # obtain the input checksums for the given file for all output files of the job
                    checksums = self.workflow.persistence.input_checksums(job, f)
                    if len(checksums) > 1:
                        # more than one checksum recorded, cannot be all the same
                        is_same = False
                    elif not checksums:
                        # no checksums recorded, we cannot assume them to be the same
                        is_same = False
                    else:
                        is_same = await f.is_same_checksum(checksums.pop())

                is_same_checksum_cache[(f, job)] = is_same
                return is_same

        def is_forced(job):
            return (
                job not in self.omitforce
                and job.rule in self.forcerules
                or not self.forcefiles.isdisjoint(job.output)
            )

        async def update_needrun(job):
            reason = self.reason(job)
            noinitreason = not reason

            if is_forced(job):
                reason.forced = True
            elif job in self.targetjobs:
                # TODO find a way to handle added/removed input files here?
                if not job.has_products(include_logfiles=False):
                    if job.input:
                        if job.rule.norun:
                            reason.updated_input_run.update(
                                [f for f in job.input if not await f.exists()]
                            )
                            if not reason.updated_input_run:
                                reason.unfinished_queue_input = (
                                    job.has_unfinished_queue_input()
                                )
                        else:
                            reason.nooutput = True
                    else:
                        reason.noio = True
                else:
                    if job.rule in self.targetrules:
                        files = set(job.products(include_logfiles=False))
                    elif (
                        self.workflow.dag_settings.target_jobs
                        and job.rule.name in self.target_jobs_rules
                    ):
                        files = set(job.products(include_logfiles=False))
                    else:
                        files = set(chain(*self.depending[job].values()))
                        if self.targetfiles:
                            files.update(
                                f for f in job.products() if f in self.targetfiles
                            )
                    reason.missing_output.update(
                        [f async for f in job.missing_output(files)]
                    )
            if not reason:
                output_mintime_ = output_mintime.get(job)
                updated_input = None
                if output_mintime_:
                    # Input is updated if it is newer that the oldest output file
                    # and does not have the same checksum as the one previously recorded.
                    updated_input = [
                        f
                        for f in job.input
                        if await f.exists()
                        and await f.is_newer(output_mintime_)
                        and not await is_same_checksum(f, job)
                    ]
                    reason.updated_input.update(updated_input)
                if not updated_input:
                    reason.unfinished_queue_input = job.has_unfinished_queue_input()
                    if not reason.unfinished_queue_input:
                        # check for other changes like parameters, set of input files, or code
                        depends_on_checkpoint_target = any(
                            f.flags.get("checkpoint_target") for f in job.input
                        )

                        if not depends_on_checkpoint_target:
                            # When the job depends on a checkpoint, it will be reevaluated in a second pass
                            # after the checkpoint output has been determined.
                            # The first pass (with depends_on_checkpoint_target == True) is not informative
                            # for determining any other changes than file modification dates, as it will
                            # change after evaluating the input function of the job in the second pass.

                            if not self.workflow.persistence.has_metadata(job):
                                reason.no_metadata = True
                            else:
                                if self.workflow.persistence.has_outdated_metadata(job):
                                    reason.outdated_metadata = True
                                if RerunTrigger.PARAMS in self.workflow.rerun_triggers:
                                    reason.params_changed = (
                                        self.workflow.persistence.params_changed(job)
                                    )
                                if RerunTrigger.INPUT in self.workflow.rerun_triggers:
                                    reason.input_changed = any(
                                        self.workflow.persistence.input_changed(job)
                                    )
                                if RerunTrigger.CODE in self.workflow.rerun_triggers:
                                    # The list comprehension is needed below in order to
                                    # collect all the async generator items before
                                    # applying any().
                                    reason.code_changed = any(
                                        [
                                            f
                                            async for f in job.outputs_older_than_script_or_notebook()
                                        ]
                                    ) or any(
                                        self.workflow.persistence.code_changed(job)
                                    )
                                if (
                                    RerunTrigger.SOFTWARE_ENV
                                    in self.workflow.rerun_triggers
                                ):
                                    reason.software_stack_changed = any(
                                        self.workflow.persistence.software_stack_changed(
                                            job
                                        )
                                    )

            if noinitreason and reason:
                reason.derived = False
            return reason

        reason = self.reason
        _needrun = self._needrun
        dependencies = self._dependencies
        depending = self.depending
        _n_until_ready = self._n_until_ready

        _needrun.clear()
        self._checkpoint_jobs.clear()
        _n_until_ready.clear()
        self._ready_jobs.clear()

        candidates = [
            [job for job in level if not self.finished(job)]
            for level in self.toposorted()
        ]

        is_all_forced = all(
            is_forced(job)
            for level in candidates
            for job in level
            if job not in self._finished
        )

        if is_all_forced:
            # all jobs are forced, no need to check for needrun
            for level in candidates:
                for job in level:
                    self.reason(job).forced = True
                    _needrun.add(job)
        else:
            # Update the output mintime of all jobs.
            # We traverse them in BFS (level order) starting from target jobs.
            # Then, we check output mintime of job itself and all direct descendants,
            # which have already been visited in the level before.
            # This way, we achieve a linear runtime.
            async with asyncio.TaskGroup() as tg:
                for level in reversed(candidates):
                    for job in level:
                        tg.create_task(update_output_mintime(job))

            # Update prior reason for all candidate jobs
            # Move from the first level to the last of the toposorted candidates.
            # If a job is needrun, mask all downstream jobs, they will below
            # in the bi-directional BFS
            # be determined as needrun because they depend on them.
            masked = set()
            queue = deque()
            for level in candidates:
                for job in level:
                    if job in masked:
                        # depending jobs of jobs that are needrun as a prior
                        # can be skipped
                        continue
                    if await update_needrun(job):
                        queue.append(job)
                        masked.update(self.bfs(self.depending, job))

            # bi-directional BFS to determine further needrun jobs
            visited = set(queue)
            candidates_set = set(job for level in candidates for job in level)
            while queue:
                job = queue.popleft()
                _needrun.add(job)

                for job_, files in dependencies[job].items():
                    if job_ in candidates_set:
                        missing_output = [f async for f in job_.missing_output(files)]
                        reason(job_).missing_output.update(missing_output)
                        if missing_output and job_ not in visited:
                            visited.add(job_)
                            queue.append(job_)

                for job_, files in depending[job].items():
                    if job_ in candidates_set:
                        if job_ not in visited:
                            if all([f.is_ancient and await f.exists() for f in files]):
                                # No other reason to run job_.
                                # Since all files are ancient, we do not trigger it.
                                continue
                            visited.add(job_)
                            queue.append(job_)
                        reason(job_).updated_input_run.update(files)

        # update _n_until_ready
        for job in _needrun:
            _n_until_ready[job] = sum(1 for dep in dependencies[job] if dep in _needrun)
            if job.is_checkpoint:
                self._checkpoint_jobs.add(job)

        # update len including finished jobs (because they have already increased the job counter)
        self._len = len(self._finished | self._needrun)

    def in_until(self, job):
        """Return whether given job has been specified via --until."""
        return job.rule.name in self.untilrules or not self.untilfiles.isdisjoint(
            job.output
        )

    def in_omitfrom(self, job):
        """Return whether given job has been specified via --omit-from."""
        return job.rule.name in self.omitrules or not self.omitfiles.isdisjoint(
            job.output
        )

    def until_jobs(self):
        """Returns a generator of jobs specified by untiljobs."""
        return (job for job in self.jobs if self.in_until(job))

    def omitfrom_jobs(self):
        """Returns a generator of jobs specified by omitfromjobs."""
        return (job for job in self.jobs if self.in_omitfrom(job))

    def downstream_of_omitfrom(self):
        """Returns the downstream of --omit-from rules or files and themselves."""
        return self.bfs(self.depending, *self.omitfrom_jobs())

    def delete_omitfrom_jobs(self):
        """Removes jobs downstream of jobs specified by --omit-from."""
        if not self.omitrules and not self.omitfiles:
            return
        downstream_jobs = list(
            self.downstream_of_omitfrom()
        )  # need to cast as list before deleting jobs
        for job in downstream_jobs:
            self.delete_job(job, recursive=False, add_dependencies=True)

    def set_until_jobs(self):
        """Removes jobs downstream of jobs specified by --omit-from."""
        if not self.untilrules and not self.untilfiles:
            return
        self.targetjobs = set(self.until_jobs())

    def update_priority(self):
        """Update job priorities."""
        prioritized = (
            lambda job: job.rule in self.priorityrules
            or not self.priorityfiles.isdisjoint(job.output)
        )
        for job in self.needrun_jobs():
            self._priority[job] = job.rule.priority
        for job in self.bfs(
            self._dependencies,
            *filter(prioritized, self.needrun_jobs()),
            stop=self.noneedrun_finished,
        ):
            self._priority[job] = Job.HIGHEST_PRIORITY

    def update_groups(self):
        groups = dict()
        for job in self.needrun_jobs():
            if job.group is None:
                continue
            stop = lambda j: j.group != job.group
            # BFS into depending needrun jobs if in same group
            # Note: never go up here (into depending), because it may contain
            # jobs that have been sorted out due to e.g. ruleorder.
            group = self.group_job_factory.new(
                job.group,
                (
                    job
                    for job in self.bfs(self._dependencies, job, stop=stop)
                    if self.needrun(job)
                ),
                self.workflow.global_resources,
            )

            # merge with previously determined groups if present
            for j in group:
                if j in groups:
                    other = groups[j]
                    other.merge(group)
                    group = other
            # update assignment
            for j in group:
                # Since groups might have been merged, we need
                # to update each job j in group.
                groups[j] = group

        self._group = groups

        self._update_group_components()

        self._check_groups()

    def _update_group_components(self):
        # span connected components if requested
        groups_by_id = defaultdict(set)
        for group in self._group.values():
            groups_by_id[group.groupid].add(group)
        for groupid, conn_components in groups_by_id.items():
            lookup_id = self.workflow.parent_groupids.get(groupid, groupid)
            n_components = self.workflow.group_settings.group_components.get(
                lookup_id, 1
            )
            if n_components > 1:
                for chunk in group_into_chunks(n_components, conn_components):
                    if len(chunk) > 1:
                        primary = chunk[0]
                        for secondary in chunk[1:]:
                            primary.merge(secondary)
                        for j in primary:
                            self._group[j] = primary

        for group in self._group.values():
            group.finalize()

    def _check_groups(self):
        """Check whether all groups are valid."""

        # find paths of jobs that leave a group and then enter it again
        # this is not allowed since then the group depends on itself
        def dfs(job, group, visited, outside_jobs, outside_jobs_all, skip_this):
            """Inner function for DFS traversal."""
            if job in group:
                if not skip_this and outside_jobs:
                    outside_jobs_all[job] = outside_jobs
                    return
            else:
                outside_jobs.append(job)
            for job_ in self._dependencies[job]:
                if job_ not in visited:
                    visited.add(job_)
                    dfs(
                        job_,
                        group,
                        visited,
                        list(outside_jobs),
                        outside_jobs_all,
                        False,
                    )

        for group in self._group.values():
            for job in group:
                outside_jobs_all = dict()
                dfs(job, group, set(), [], outside_jobs_all, True)
                if outside_jobs_all:
                    fmt_outside = lambda jobs: ",".join(
                        sorted(set(j.rule.name for j in jobs))
                    )
                    bullet = "* " if len(outside_jobs_all) > 1 else ""
                    fixes = "\n".join(
                        f"{bullet}Remove {job.rule.name} from the group or add {fmt_outside(outside)} to the group."
                        for job, outside in outside_jobs_all.items()
                    )
                    raise WorkflowError(
                        f"Group {group.groupid} depends on itself. "
                        "This is not allowed, because it would lead to an "
                        "the group can never be ready for execution. "
                        "Ensure that there is no path of jobs in the DAG that "
                        "starts in a group, leaves it (i.e. at least one job in "
                        "the path is not in the group), and then enters it again. "
                        f"Possible fixes are:\n{fixes}",
                        rule=job.rule,
                    )

    async def update_incomplete_input_expand_jobs(self):
        """Update (re-evaluate) all jobs which have incomplete input file expansions.

        only filled in the second pass of postprocessing.
        """
        updated = False
        for job in list(self.jobs):
            if job.incomplete_input_expand:
                newjob = await job.updated()
                await self.replace_job(job, newjob, recursive=False)
                updated = True
        return updated

    def update_ready(self, jobs=None):
        """Update information whether a job is ready to execute.

        Given jobs must be needrun jobs!
        """

        if jobs is None:
            jobs = self.needrun_jobs()

        potential_new_ready_jobs = False
        candidate_groups = set()
        for job in jobs:
            if job in self._ready_jobs or job in self._running:
                # job has been seen before or is running, no need to process again
                continue
            if (
                not self.finished(job)
                and self._ready(job)
                and not job.has_unfinished_queue_input()
            ):
                potential_new_ready_jobs = True
                if job.group is None:
                    self._ready_jobs.add(job)
                else:
                    group = self._group[job]
                    if group not in self._running:
                        candidate_groups.add(group)

        self._ready_jobs.update(
            group
            for group in candidate_groups
            if all(self._ready(job) for job in group)
        )
        return potential_new_ready_jobs

    def get_jobs_or_groups(self):
        visited_groups = set()
        for job in self.jobs:
            if job.group is None:
                yield job
            else:
                group = self._group[job]
                if group in visited_groups:
                    continue
                visited_groups.add(group)
                yield group

    async def postprocess(
        self, update_needrun=True, update_incomplete_input_expand_jobs=True
    ):
        """Postprocess the DAG. This has to be invoked after any change to the
        DAG topology."""
        self.cleanup()
        self.update_jobids()
        if update_needrun:
            # cleanup local storage copies of files that will be created by jobs
            # this is important to ensure that there are no outdated local copies
            # that misguide e.g. params functions.
            await self.sanitize_local_storage_copies()
            self.update_container_imgs()
            self.update_conda_envs()
            await self.update_needrun()
        self.update_priority()
        self.handle_pipes_and_services()
        self.handle_update_flags()
        self.update_groups()
        self.update_storage_inputs()

        if update_incomplete_input_expand_jobs:
            updated = await self.update_incomplete_input_expand_jobs()
            if updated:
                # run a second pass, some jobs have been updated
                # with potentially new input files that have depended
                # on group ids.
                await self.postprocess(
                    update_needrun=True, update_incomplete_input_expand_jobs=False
                )

                return

        self.update_ready()
        self.update_checkpoint_outputs()

    def handle_pipes_and_services(self):
        """Use pipes and services to determine job groups. Check if every pipe has exactly
        one consumer"""

        visited = set()
        for job in self.needrun_jobs():
            candidate_groups = set()
            user_groups = set()
            if job.pipe_group is not None:
                candidate_groups.add(job.pipe_group)
            if job.group is not None:
                user_groups.add(job.group)
            all_depending = set()
            has_pipe_or_service = False
            for f in job.output:
                is_pipe = is_flagged(f, "pipe")
                is_service = is_flagged(f, "service")
                if is_pipe or is_service:
                    if job.is_run:
                        raise WorkflowError(
                            "Rule defines pipe output but "
                            "uses a 'run' directive. This is "
                            "not possible for technical "
                            "reasons. Consider using 'shell' or "
                            "'script'.",
                            rule=job.rule,
                        )

                    has_pipe_or_service = True
                    depending = [
                        j for j, files in self.depending[job].items() if f in files
                    ]
                    if is_pipe and len(depending) > 1:
                        raise WorkflowError(
                            "Output file {} is marked as pipe "
                            "but more than one job depends on "
                            "it. Make sure that any pipe "
                            "output is only consumed by one "
                            "job".format(f),
                            rule=job.rule,
                        )
                    elif len(depending) == 0:
                        raise WorkflowError(
                            "Output file {} is marked as pipe or service "
                            "but it has no consumer. This is "
                            "invalid because it can lead to "
                            "a dead lock.".format(f),
                            rule=job.rule,
                        )
                    elif is_pipe and depending[0].is_norun:
                        raise WorkflowError(
                            f"Output file {f} is marked as pipe but is requested by a rule that "
                            "does not execute anything. This is not allowed because it would lead "
                            "to a dead lock."
                        )

                    for dep in depending:
                        if dep.is_run:
                            raise WorkflowError(
                                "Rule consumes pipe or service input but "
                                "uses a 'run' directive. This is "
                                "not possible for technical "
                                "reasons. Consider using 'shell' or "
                                "'script'.",
                                rule=job.rule,
                            )

                        all_depending.add(dep)
                        if dep.pipe_group is not None:
                            candidate_groups.add(dep.pipe_group)
                        if dep.group is not None:
                            user_groups.add(dep.group)

            if not has_pipe_or_service:
                continue

            # All pipe groups should be contained within one user-defined group
            if len(user_groups) > 1:
                raise WorkflowError(
                    "An output file is marked as "
                    "pipe or service, but consuming jobs "
                    "are part of conflicting "
                    "groups.",
                    rule=job.rule,
                )

            if len(candidate_groups) > 1:
                # Merge multiple pipe groups together
                group = candidate_groups.pop()
                for g in candidate_groups:
                    g.merge(group)
            elif candidate_groups:
                # extend the candidate group to all involved jobs
                group = candidate_groups.pop()
            else:
                # generate a random unique group name
                group = CandidateGroup()  # str(uuid.uuid4())

            # Assign the pipe group to all involved jobs.
            job.pipe_group = group
            visited.add(job)
            for j in all_depending:
                j.pipe_group = group
                visited.add(j)

        # convert candidate groups to plain string IDs
        for job in visited:
            # Set the group every job with an assigned pipe_group but no user-defined
            # group to the pipe_group
            if job.pipe_group and job.group is None:
                job.group = job.pipe_group.id

    def _ready(self, job):
        """Return whether the given job is ready to execute."""
        group = self._group.get(job, None)
        if group is None:
            return self._n_until_ready[job] == 0
        else:
            n_internal_deps = lambda job: sum(
                self._group.get(dep) == group for dep in self._dependencies[job]
            )
            return all(
                (self._n_until_ready[job] - n_internal_deps(job)) == 0 for job in group
            )

    async def update_queue_input_jobs(self):
        updated = False
        if self.has_unfinished_queue_input_jobs():
            logger.info("Updating jobs with queue input...")
            for job in self.queue_input_jobs:
                if (
                    job.has_queue_input()
                    and job not in self._jobs_with_finished_queue_input
                ):
                    newjob = await job.updated()
                    if newjob.input != job.input:
                        await self.replace_job(job, newjob, recursive=False)
                        updated = True
                    if updated and not job.has_unfinished_queue_input():
                        self._jobs_with_finished_queue_input.add(job)
            if updated:
                await self.postprocess_after_update()
                # reset queue_input_jobs such that it is recomputed next time
                self._queue_input_jobs = None
        return updated

    @property
    def queue_input_jobs(self):
        if self._queue_input_jobs is None:
            self._queue_input_jobs = set(
                job
                for job in self.needrun_jobs(exclude_finished=False)
                if job.has_queue_input()
            )
        return self._queue_input_jobs

    def has_unfinished_queue_input_jobs(self):
        return any(job.has_unfinished_queue_input() for job in self.queue_input_jobs)

    async def update_checkpoint_dependencies(self, jobs=None):
        """Update dependencies of checkpoints."""
        updated = False
        self.update_checkpoint_outputs()
        if jobs is None:
            jobs = [job for job in self.jobs if not self.needrun(job)]
        all_depending = []
        for job in jobs:
            if job.is_checkpoint:
                depending = list(self.depending[job])
                all_depending.extend(depending)
        for j in all_depending:
            newjob = await j.updated()
            await self.replace_job(j, newjob, recursive=False)
            updated = True
        if updated:
            await self.postprocess_after_update()
        return updated

    async def postprocess_after_update(self):
        await self.postprocess()
        shared_input_output = (
            SharedFSUsage.INPUT_OUTPUT in self.workflow.storage_settings.shared_fs_usage
        )
        if not self.workflow.dryrun and (
            (self.workflow.is_main_process and shared_input_output)
            or self.workflow.remote_exec
        ):
            await self.retrieve_storage_inputs()

    def register_running(self, jobs):
        self._running.update(jobs)
        self._ready_jobs -= jobs
        for job in jobs:
            try:
                del self._n_until_ready[job]
            except KeyError:
                # already gone
                pass

    async def finish(self, job, update_checkpoint_dependencies=True):
        """Finish a given job (e.g. remove from ready jobs, mark depending jobs
        as ready)."""

        self._running.remove(job)

        # turn off this job's Reason
        if job.is_group():
            for j in job:
                self.reason(j).mark_finished()
        else:
            self.reason(job).mark_finished()

        try:
            self._ready_jobs.remove(job)
        except KeyError:
            pass

        if job.is_group():
            jobs = job
        else:
            jobs = [job]

        self._finished.update(jobs)
        self.checkpoint_jobs.difference_update(job for job in jobs if job.is_checkpoint)

        updated_dag = False
        if update_checkpoint_dependencies:
            updated_dag = await self.update_checkpoint_dependencies(jobs)

        depending = [
            j
            for job in jobs
            for j in self.depending[job]
            if not self.in_until(job) and self.needrun(j)
        ]

        if not updated_dag:
            # Mark depending jobs as ready.
            # Skip jobs that are marked as until jobs.
            # This is not necessary if the DAG has been fully updated above.
            for job in depending:
                self._n_until_ready[job] -= 1

        potential_new_ready_jobs = self.update_ready(depending)

        if updated_dag:
            # We might have new jobs, so we need to ensure that all conda envs
            # and singularity images are set up.
            if (
                DeploymentMethod.APPTAINER
                in self.workflow.deployment_settings.deployment_method
            ):
                self.pull_container_imgs()
            if (
                DeploymentMethod.CONDA
                in self.workflow.deployment_settings.deployment_method
            ):
                self.create_conda_envs()
            potential_new_ready_jobs = True

        if not self.checkpoint_jobs:
            # While there are still checkpoint jobs, we cannot safely delete
            # temp files.
            # TODO: we maybe could be more accurate and determine whether there is a
            # checkpoint that depends on the temp file.
            for job in jobs:
                await self.handle_temp(job)

        return potential_new_ready_jobs

    async def new_job(
        self, rule, targetfile=None, format_wildcards=None, wildcards_dict=None
    ):
        """Create new job for given rule and (optional) targetfile.
        This will reuse existing jobs with the same wildcards."""
        product = rule.get_some_product()
        if targetfile is None and wildcards_dict is not None and product is not None:
            # no targetfile given, but wildcards_dict is given, hence this job seems to contain wildcards
            # just take one targetfile
            try:
                targetfile = product.apply_wildcards(wildcards_dict)
            except WildcardError as e:
                raise WorkflowError(
                    f"Given wildcards for rule {rule.name} do not match output file {repr(rule.output[0])}: {e}"
                )

        key = (rule, targetfile)
        if key in self.job_cache:
            assert targetfile is not None
            return self.job_cache[key]
        wildcards_dict = rule.get_wildcards(targetfile, wildcards_dict=wildcards_dict)
        job = await self.job_factory.new(
            rule,
            self,
            wildcards_dict=wildcards_dict,
            format_wildcards=format_wildcards,
            targetfile=targetfile,
        )
        self.cache_job(job)
        return job

    def cache_job(self, job):
        for f in job.products():
            self.job_cache[(job.rule, f)] = job

    def delete_job(self, job, recursive=True, add_dependencies=False):
        """Delete given job from DAG."""
        if job in self.targetjobs:
            self.targetjobs.remove(job)
        if add_dependencies:
            for _job in self._dependencies[job]:
                self.targetjobs.add(_job)
        for job_ in self.depending[job]:
            del self._dependencies[job_][job]
        del self.depending[job]
        for job_ in self._dependencies[job]:
            depending = self.depending[job_]
            del depending[job]
            if not depending and recursive:
                self.delete_job(job_)
        del self._dependencies[job]
        if job in self._reason:
            del self._reason[job]
        if job in self._needrun:
            self._len -= 1
            self._needrun.remove(job)
        if job in self._finished:
            self._finished.remove(job)
        if job in self._ready_jobs:
            self._ready_jobs.remove(job)
        if job in self._n_until_ready:
            del self._n_until_ready[job]
        # remove from cache
        for f in job.output:
            try:
                del self.job_cache[(job.rule, f)]
            except KeyError:
                pass

    async def replace_job(self, job, newjob, recursive=True):
        """Replace given job with new job."""
        add_to_targetjobs = job in self.targetjobs
        try:
            jobid = self.jobid(job)
        except KeyError:
            # Job has been added while updating another checkpoint,
            # jobid is not yet known.
            jobid = None

        depending = list(self.depending[job].items())
        if self.finished(job):
            self._finished.add(newjob)

        self.delete_job(job, recursive=recursive)
        if jobid is not None:
            self._jobid[newjob] = jobid

        if add_to_targetjobs:
            self.targetjobs.add(newjob)

        self.cache_job(newjob)

        await self.update([newjob])

        for job_, files in depending:
            self._dependencies[job_][newjob].update(files)
            self.depending[newjob][job_].update(files)

    def specialize_rule(self, rule, newrule):
        """Specialize the given rule by inserting newrule into the DAG."""
        assert newrule is not None
        self.rules.add(newrule)
        self.update_output_index()

    def is_batch_rule(self, rule):
        """Return True if the underlying rule is to be used for batching the DAG."""
        return self.batch is not None and rule.name == self.batch.rulename

    def handle_update_flags(self):
        before_update_jobs = dict()
        update_jobs = dict()
        for job in self.needrun_jobs():
            for f in job.input:
                if is_flagged(f, "before_update"):
                    before_update_jobs[f] = job
            for f in job.output:
                if is_flagged(f, "update"):
                    update_jobs[f] = job

        for f, job in before_update_jobs.items():
            update_job = update_jobs.get(f)
            if update_job is not None and job.priority <= update_job.priority:
                logger.info(
                    f"Raising priority of job {job} such that it runs before "
                    f"job {update_job} because of flag 'before_update' on {f}"
                )
                self._priority[job] = update_job.priority + 1
                f_hash = hashlib.sha256()
                f_hash.update(f.encode())
                mutex = f"update_file_{f_hash.hexdigest()}"
                job.add_aux_resource(mutex, 1)
                update_job.add_aux_resource(mutex, 1)
                self.workflow.register_resource(mutex, 1)

    async def collect_potential_dependencies(self, job, known_producers):
        """Collect all potential dependencies of a job. These might contain
        ambiguities. The keys of the returned dict represent the files to be considered.
        """
        # use a set to circumvent multiple jobs for the same file
        # if user specified it twice
        file2jobs = self.file2jobs

        input_files = list(job.unique_input)

        if self.is_batch_rule(job.rule):
            # only consider the defined partition of the input files
            input_batch = self.batch.get_batch(input_files)
            if len(input_batch) != len(input_files):
                logger.info(
                    "Considering only batch {} for DAG computation.\n"
                    "All jobs beyond the batching rule are omitted until the final batch.\n"
                    "Don't forget to run the other batches too.".format(self.batch)
                )
                input_files = input_batch

        for file in input_files:
            if is_flagged(file, "before_update"):
                # do not find a producer for this file, it shall be considered in its
                # form before the update
                yield PotentialDependency(file, None, False)
            else:
                try:
                    yield PotentialDependency(file, known_producers[file], True)
                except KeyError:
                    try:
                        if file in job.dependencies:
                            yield PotentialDependency(
                                file,
                                [
                                    await self.new_job(
                                        job.dependencies[file],
                                        targetfile=file,
                                        wildcards_dict=job.wildcards_dict,
                                    )
                                ],
                                False,
                            )
                        else:
                            yield PotentialDependency(
                                file,
                                await file2jobs(
                                    file, wildcards_dict=job.wildcards_dict
                                ),
                                False,
                            )
                    except MissingRuleException as ex:
                        # no dependency found
                        yield PotentialDependency(file, None, False)

    def bfs(self, direction, *jobs, stop=lambda job: False):
        """Perform a breadth-first traversal of the DAG."""
        queue = deque(jobs)
        visited = set(queue)
        while queue:
            job = queue.popleft()
            if stop(job):
                # stop criterion reached for this node
                continue
            yield job
            for job_ in direction[job].keys():
                if not job_ in visited:
                    queue.append(job_)
                    visited.add(job_)

    def level_bfs(self, direction, *jobs, stop=lambda job: False):
        """Perform a breadth-first traversal of the DAG, but also yield the
        level together with each job."""
        queue = [(job, 0) for job in jobs]
        visited = set(jobs)
        while queue:
            job, level = queue.pop(0)
            if stop(job):
                # stop criterion reached for this node
                continue
            yield level, job
            level += 1
            for job_, _ in direction[job].items():
                if not job_ in visited:
                    queue.append((job_, level))
                    visited.add(job_)

    def dfs(self, direction, *jobs, stop=lambda job: False, post=True):
        """Perform depth-first traversal of the DAG."""
        visited = set()

        def _dfs(job):
            """Inner function for DFS traversal."""
            if stop(job):
                return
            if not post:
                yield job
            for job_ in direction[job]:
                if job_ not in visited:
                    visited.add(job_)
                    for j in _dfs(job_):
                        yield j
            if post:
                yield job

        for job in jobs:
            for job_ in _dfs(job):
                yield job_

    def new_wildcards(self, job):
        """Return wildcards that are newly introduced in this job,
        compared to its ancestors."""
        new_wildcards = set(job.wildcards.items())
        for job_ in self._dependencies[job]:
            if not new_wildcards:
                return set()
            for wildcard in job_.wildcards.items():
                new_wildcards.discard(wildcard)
        return new_wildcards

    async def rule2job(self, targetrule):
        """Generate a new job from a given rule."""
        if targetrule.has_wildcards():
            raise WorkflowError(
                "Target rules may not contain wildcards. "
                "Please specify concrete files or a rule without wildcards at the command line, "
                "or have a rule without wildcards at the very top of your workflow (e.g. the typical "
                '"rule all" which just collects all results you want to generate in the end).'
            )
        return await self.new_job(targetrule)

    async def file2jobs(self, targetfile, wildcards_dict=None):
        rules = self.output_index.match(targetfile)
        jobs = []
        exceptions = list()
        for rule in rules:
            if rule.is_producer(targetfile):
                try:
                    jobs.append(
                        await self.new_job(
                            rule, targetfile=targetfile, wildcards_dict=wildcards_dict
                        )
                    )
                except InputFunctionException as e:
                    exceptions.append(e)
        if not jobs:
            if exceptions:
                raise exceptions[0]
            raise MissingRuleException(targetfile)
        return jobs

    def rule_dot2(self):
        dag = defaultdict(list)
        visited = set()
        preselect = set()

        def preselect_parents(job):
            for parent in self.depending[job]:
                if parent in preselect:
                    continue
                preselect.add(parent)
                preselect_parents(parent)

        def build_ruledag(job, key=lambda job: job.rule.name):
            if job in visited:
                return
            visited.add(job)
            deps = sorted(self._dependencies[job], key=key)
            deps = [
                (
                    group[0]
                    if preselect.isdisjoint(group)
                    else preselect.intersection(group).pop()
                )
                for group in (list(g) for _, g in groupby(deps, key))
            ]
            dag[job].extend(deps)
            preselect_parents(job)
            for dep in deps:
                build_ruledag(dep)

        for job in self.targetjobs:
            build_ruledag(job)

        return self._dot(dag.keys())

    def rule_dot(self):
        graph = defaultdict(set)
        for job in self.jobs:
            graph[job.rule].update(dep.rule for dep in self._dependencies[job])
        return self._dot(graph)

    def dot(self):
        def node2style(job):
            if not self.needrun(job):
                return "rounded,dashed"
            return "rounded"

        def format_wildcard(wildcard):
            name, value = wildcard
            return f"{name}: {value}"

        node2rule = lambda job: job.rule
        node2label = lambda job: "\\n".join(
            chain(
                [job.rule.name], sorted(map(format_wildcard, self.new_wildcards(job)))
            )
        )

        dag = {job: self._dependencies[job] for job in self.jobs}

        return self._dot(
            dag, node2rule=node2rule, node2style=node2style, node2label=node2label
        )

    def _dot(
        self,
        graph,
        node2rule=lambda node: node,
        node2style=lambda node: "rounded",
        node2label=lambda node: node,
    ):
        # color rules
        huefactor = 2 / (3 * len(self.rules))
        rulecolor = {
            rule: "{:.2f} 0.6 0.85".format(i * huefactor)
            for i, rule in enumerate(self.rules)
        }

        # markup
        node_markup = '\t{}[label = "{}", color = "{}", style="{}"];'.format
        edge_markup = "\t{} -> {}".format

        # node ids
        ids = {node: i for i, node in enumerate(graph)}

        # calculate nodes
        nodes = [
            node_markup(
                ids[node],
                node2label(node),
                rulecolor[node2rule(node)],
                node2style(node),
            )
            for node in graph
        ]
        # calculate edges
        edges = [
            edge_markup(ids[dep], ids[node])
            for node, deps in graph.items()
            for dep in deps
        ]

        return textwrap.dedent(
            """\
            digraph snakemake_dag {{
                graph[bgcolor=white, margin=0];
                node[shape=box, style=rounded, fontname=sans, \
                fontsize=10, penwidth=2];
                edge[penwidth=2, color=grey];
            {items}
            }}\
            """
        ).format(items="\n".join(nodes + edges))

    def filegraph_dot(
        self,
        node2rule=lambda node: node,
        node2style=lambda node: "rounded",
        node2label=lambda node: node,
    ):
        # NOTE: This is code from the rule_dot method.
        # This method could be split like there as well, however,
        # it cannot easily reuse the _dot method due to the different node type
        graph = defaultdict(set)
        for job in self.jobs:
            graph[job.rule].update(dep.rule for dep in self._dependencies[job])

        # node ids
        ids = {node: i for i, node in enumerate(graph)}

        # Compute colors for rules
        def hsv_to_htmlhexrgb(h, s, v):
            """Convert hsv colors to hex-encoded rgb colors usable by html."""
            import colorsys

            hex_r, hex_g, hex_b = (round(255 * x) for x in colorsys.hsv_to_rgb(h, s, v))
            return "#{hex_r:0>2X}{hex_g:0>2X}{hex_b:0>2X}".format(
                hex_r=hex_r, hex_g=hex_g, hex_b=hex_b
            )

        huefactor = 2 / (3 * len(self.rules))
        rulecolor = {
            rule: hsv_to_htmlhexrgb(i * huefactor, 0.6, 0.85)
            for i, rule in enumerate(self.rules)
        }

        def resolve_input_functions(input_files):
            """Iterate over all input files and replace input functions
            with a fixed string.
            """
            files = []
            for f in input_files:
                if callable(f):
                    files.append("<input function>")
                    # NOTE: This is a workaround. It would be more informative
                    # to show the code of the input function here (if it is
                    # short enough). This cannot be easily done with the inspect
                    # module, since the line numbers in the Snakefile do not
                    # behave as expected. One (complicated) solution for this
                    # would be to find the Snakefile and directly extract the
                    # code of the function.
                else:
                    files.append(repr(f).strip("'"))
            return files

        def html_node(node_id, node, color):
            """Assemble a html style node for graphviz"""
            input_files = resolve_input_functions(node._input)
            output_files = [repr(f).strip("'") for f in node._output]
            input_header = (
                '<b><font point-size="14">&#8618; input</font></b>'
                if input_files
                else ""
            )
            output_header = (
                '<b><font point-size="14">output &rarr;</font></b>'
                if output_files
                else ""
            )
            html_node = [
                '{node_id} [ shape=none, margin=0, label=<<table border="2" color="{color}" cellspacing="3" cellborder="0">'.format(
                    node_id=node_id, color=color
                ),
                "<tr><td>",
                f'<b><font point-size="18">{node.name}</font></b>',
                "</td></tr>",
                "<hr/>",
                '<tr><td align="left"> {input_header} </td></tr>'.format(
                    input_header=input_header
                ),
            ]

            for filename in sorted(input_files):
                # Escape html relevant chars like '<' and '>' in filenames
                # These can be added by input functions etc. and cannot be
                # displayed in graphviz HTML nodes.
                in_file = html.escape(filename)
                html_node.extend(
                    [
                        "<tr>",
                        '<td align="left"><font face="monospace">{in_file}</font></td>'.format(
                            in_file=in_file
                        ),
                        "</tr>",
                    ]
                )

            html_node.append("<hr/>")
            html_node.append(
                '<tr><td align="right"> {output_header} </td> </tr>'.format(
                    output_header=output_header
                )
            )

            for filename in sorted(output_files):
                out_file = html.escape(filename)
                html_node.extend(
                    [
                        "<tr>",
                        '<td align="left"><font face="monospace">{out_file}</font></td>'
                        "</tr>".format(out_file=out_file),
                    ]
                )

            html_node.append("</table>>]")
            return "\n".join(html_node)

        nodes = [
            html_node(ids[node], node, rulecolor[node2rule(node)]) for node in graph
        ]

        # calculate edges
        edge_markup = "\t{} -> {}".format
        edges = [
            edge_markup(ids[dep], ids[node], ids[dep], ids[node])
            for node, deps in graph.items()
            for dep in deps
        ]

        return textwrap.dedent(
            """\
            digraph snakemake_dag {{
                graph[bgcolor=white, margin=0];
                node[shape=box, style=rounded, fontname=sans, \
                fontsize=10, penwidth=2];
                edge[penwidth=2, color=grey];
            {items}
            }}\
            """
        ).format(items="\n".join(nodes + edges))

    async def summary(self, detailed=False):
        def fmt_output(f):
            if f.is_storage:
                return f.storage_object.query
            return f

        if detailed:
            yield "output_file\tdate\trule\tlog-file(s)\tinput-file(s)\tshellcmd\tstatus\tplan"
        else:
            yield "output_file\tdate\trule\tlog-file(s)\tstatus\tplan"

        for job in self.jobs:
            output = job.output
            for f in output:
                rule = job.rule.name

                date = (
                    time.ctime((await f.mtime()).local_or_storage())
                    if await f.exists()
                    else "-"
                )

                pending = "update pending" if self.reason(job) else "no update"

                log = self.workflow.persistence.log(f)
                log = "-" if log is None else ",".join(log)

                inputfiles = self.workflow.persistence.input(f)
                inputfiles = "-" if inputfiles is None else ",".join(inputfiles)

                shellcmd = self.workflow.persistence.shellcmd(f)
                shellcmd = "-" if shellcmd is None else shellcmd
                # remove new line characters, leading and trailing whitespace
                shellcmd = shellcmd.strip().replace("\n", "; ")

                status = "ok"
                if not await f.exists():
                    if is_flagged(f, "temp"):
                        status = "removed temp file"
                    elif is_flagged(f, "pipe"):
                        status = "pipe file"
                    elif is_flagged(f, "service"):
                        status = "service file"
                    else:
                        status = "missing"
                elif self.reason(job).updated_input:
                    status = "updated input files"
                elif self.workflow.persistence.code_changed(job, file=f):
                    status = "rule implementation changed"
                elif self.workflow.persistence.input_changed(job, file=f):
                    status = "set of input files changed"
                elif self.workflow.persistence.params_changed(job, file=f):
                    status = "params changed"
                if detailed:
                    yield "\t".join(
                        (
                            fmt_output(f),
                            date,
                            rule,
                            log,
                            inputfiles,
                            shellcmd,
                            status,
                            pending,
                        )
                    )
                else:
                    yield "\t".join((fmt_output(f), date, rule, log, status, pending))

    def archive(self, path: Path):
        """Archives workflow such that it can be re-run on a different system.

        Archiving includes git versioned files (i.e. Snakefiles, config files, ...),
        ancestral input files and conda environments.
        """
        if path.suffix == ".tar":
            mode = "x"
        elif path.suffixes == [".tar", ".bz2"]:
            mode = "x:bz2"
        elif path.suffixes == [".tar", ".xz"]:
            mode = "x:xz"
        elif path.suffixes == [".tar", ".gz"]:
            mode = "x:gz"
        else:
            raise WorkflowError(
                "Unsupported archive format "
                "(supported: .tar, .tar.gz, .tar.bz2, .tar.xz)"
            )
        if path.exists():
            raise WorkflowError(f"Archive already exists:\n{path}")

        self.create_conda_envs()

        try:
            workdir = Path(os.path.abspath(os.getcwd()))
            with tarfile.open(path, mode=mode, dereference=True) as archive:
                archived = set()

                def add(path):
                    if workdir not in Path(os.path.abspath(path)).parents:
                        logger.warning(
                            "Path {} cannot be archived: "
                            "not within working directory.".format(path)
                        )
                    else:
                        f = os.path.relpath(path)
                        if f not in archived:
                            archive.add(f)
                            archived.add(f)
                            logger.info("archived " + f)

                logger.info(
                    "Archiving snakefiles, scripts and files under "
                    "version control..."
                )
                for f in self.get_sources():
                    add(f)

                logger.info("Archiving external input files...")
                for job in self.jobs:
                    # input files
                    for f in job.input:
                        if self.is_external_input(f, job):
                            # this is an input file that is not created by any job
                            add(f)

                logger.info("Archiving conda environments...")
                envs = set()
                for job in self.jobs:
                    if job.conda_env_spec:
                        env_archive = job.archive_conda_env()
                        envs.add(env_archive)
                for env in envs:
                    add(env)

        except BaseException as e:
            os.remove(path)
            raise e

    def is_external_input(self, file, job, not_needrun_is_external=False):
        """Return True if the given file is an external input for the given job."""
        consider = lambda job: True
        if not_needrun_is_external:
            consider = lambda job: self.needrun(job)
        return not any(
            file in files
            for dep, files in self._dependencies[job].items()
            if consider(dep)
        )

    async def clean(self, only_temp=False, dryrun=False):
        """Removes files generated by the workflow."""
        for job in self.jobs:
            for f in job.output:
                if not only_temp or is_flagged(f, "temp"):
                    # The reason for the second check is that dangling
                    # symlinks fail f.exists.
                    if await f.exists() or os.path.islink(f):
                        if await f.protected():
                            logger.error(f"Skipping write-protected file {f}.")
                        else:
                            msg = "Deleting {}" if not dryrun else "Would delete {}"
                            logger.info(msg.format(f))
                            if not dryrun:
                                # Remove non-empty dirs if flagged as temp()
                                await f.remove(remove_non_empty_dir=only_temp)

    def list_untracked(self):
        """List files in the workdir that are not in the dag."""
        used_files = set()
        files_in_cwd = set()
        for job in self.jobs:
            used_files.update(
                os.path.relpath(file)
                for file in chain(job.local_input, job.local_output, job.log)
            )
        for root, dirs, files in os.walk(os.getcwd()):
            # Ignore hidden files and don't traverse into hidden dirs
            files_in_cwd.update(
                [
                    os.path.relpath(os.path.join(root, f))
                    for f in files
                    if not f[0] == "."
                ]
            )
            dirs[:] = [d for d in dirs if not d[0] == "."]
        for f in sorted(list(files_in_cwd - used_files)):
            print(f)

    def d3dag(self, max_jobs=10000):
        def node(job):
            jobid = self.jobid(job)
            return {
                "id": jobid,
                "value": {
                    "jobid": jobid,
                    "label": job.rule.name,
                    "rule": job.rule.name,
                },
            }

        def edge(a, b):
            return {"u": self.jobid(a), "v": self.jobid(b)}

        jobs = list(self.jobs)

        if len(jobs) > max_jobs:
            logger.info(f"Job-DAG is too large for visualization (>{max_jobs} jobs).")
        else:
            logger.d3dag(
                nodes=[node(job) for job in jobs],
                edges=[
                    edge(dep, job)
                    for job in jobs
                    for dep in self._dependencies[job]
                    if self.needrun(dep)
                ],
            )

    def print_reasons(self):
        """Print summary of execution reasons."""
        reasons = defaultdict(set)
        for job in self.needrun_jobs(exclude_finished=False):
            for reason in self.reason(job).get_names():
                reasons[reason].add(job.rule.name)
        if reasons:
            msg = "Reasons:\n    (check individual jobs above for details)"
            for reason, rules in sorted(reasons.items()):
                rules = sorted(rules)
                if len(rules) > 50:
                    rules = rules[:50] + ["..."]
                rules = ", ".join(rules)
                msg += f"\n    {reason}:\n        {rules}"
            logger.info(msg)

    def stats(self):
        from tabulate import tabulate

        rules = Counter()
        rules.update(job.rule for job in self.needrun_jobs())
        rules.update(job.rule for job in self.finished_jobs)

        rows = [
            {
                "job": rule.name,
                "count": count,
            }
            for rule, count in sorted(
                rules.most_common(), key=lambda item: item[0].name
            )
        ]
        rows.append(
            {
                "job": "total",
                "count": sum(rules.values()),
            }
        )

        yield "Job stats:"
        yield tabulate(rows, headers="keys")
        yield ""

    def toposorted(self, jobs=None, inherit_pipe_dependencies=False):
        if jobs is None:
            jobs = set(self.jobs)

        pipe_dependencies = {}
        # If enabled, toposort should put all pipe jobs at the same level. We do this by
        # causing all jobs in the pipe group to use each other's dependencies
        if inherit_pipe_dependencies:
            # First, we organize all the jobs in the group into a dict according to
            # their pipe_group
            pipe_groups = defaultdict(list)
            for name, group in groupby(jobs, attrgetter("pipe_group")):
                if name is not None:
                    pipe_groups[name].extend(group)

            # Then, for each pipe_group, we find the dependencies of every job in the
            # group, filtering out any dependencies that are, themselves, in the group
            for name, group in pipe_groups.items():
                pipe_dependencies[name] = set(
                    d
                    for job in group
                    for d in self._dependencies[job]
                    if d not in group
                )

        # Collect every job's dependencies into a definitive mapping
        dependencies = {}
        for job in jobs:
            if job.pipe_group in pipe_dependencies:
                deps = pipe_dependencies[job.pipe_group]
            else:
                deps = self._dependencies[job]
            dependencies[job] = {dep for dep in deps if dep in jobs}

        toposorted = toposort(dependencies)

        # Within each toposort layer, entries should be sorted so that pipe jobs are
        # listed order of dependence, i.e. dependent jobs before depending jobs
        for layer in toposorted:
            pipe_groups = defaultdict(set)
            sorted_layer = []
            for job in layer:
                if job.pipe_group is None:
                    sorted_layer.append(job)
                    continue
                pipe_groups[job.pipe_group].add(job)

            for group in pipe_groups.values():
                sorted_layer.extend(
                    chain.from_iterable(
                        toposort(
                            {
                                job: {
                                    dep
                                    for dep in self._dependencies[job]
                                    if dep in group
                                }
                                for job in group
                            }
                        )
                    )
                )

            yield sorted_layer

    async def get_outputs_with_changes(
        self, change_type: ChangeType, include_needrun=True
    ):
        is_changed = lambda job: (
            getattr(self.workflow.persistence, f"{change_type}_changed")(job)
            if not job.is_group() and (include_needrun or not self.needrun(job))
            else []
        )
        changed = list(chain(*map(is_changed, self.jobs)))
        if change_type == ChangeType.CODE:
            for job in self.jobs:
                if not job.is_group() and (include_needrun or not self.needrun(job)):
                    changed.extend(
                        list(await job.outputs_older_than_script_or_notebook())
                    )
        return changed

    def get_sources(self):
        files = set()

        def local_path(f):
            if not isinstance(f, SourceFile) and is_local_file(f):
                return f
            if isinstance(f, LocalSourceFile):
                return f.get_path_or_uri()

        def norm_rule_relpath(f, rule):
            if not os.path.isabs(f):
                f = os.path.join(rule.basedir, f)
            return os.path.relpath(f)

        # get registered sources
        for f in self.workflow.included:
            f = local_path(f)
            if f:
                try:
                    f = os.path.relpath(f)
                except ValueError:
                    if ON_WINDOWS:
                        pass  # relpath doesn't work on win if files are on different drive
                    else:
                        raise
                files.add(f)
        for rule in self.workflow.rules:
            script_path = rule.script or rule.notebook
            if script_path:
                script_path = norm_rule_relpath(script_path, rule)
                files.add(script_path)
                script_dir = os.path.dirname(script_path)
                files.update(
                    os.path.join(dirpath, f)
                    for dirpath, _, files in os.walk(script_dir)
                    for f in files
                )

        for job in self.jobs:
            assert not job.is_group(), "bug: groups should not be yielded by DAG.jobs"
            if job.conda_env_spec and job.conda_env_spec.is_file:
                f = local_path(job.conda_env_spec.file)
                if f:
                    # url points to a local env file
                    env_path = norm_rule_relpath(f, job.rule)
                    files.add(env_path)

        for f in self.workflow.configfiles:
            files.add(f)

        # get git-managed files
        # TODO allow a manifest file as alternative
        try:
            out = subprocess.check_output(
                ["git", "ls-files", "--recurse-submodules", "."], stderr=subprocess.PIPE
            )
            for f in out.decode().split("\n"):
                if f:
                    files.add(os.path.relpath(f))
        except subprocess.CalledProcessError as e:
            if "fatal: not a git repository" in e.stderr.decode().lower():
                logger.warning(
                    "Unable to retrieve additional files from git. "
                    "This is not a git repository."
                )
            else:
                raise WorkflowError(
                    "Error executing git (Snakemake requires git to be installed for "
                    "remote execution without shared filesystem):\n" + e.stderr.decode()
                )

        return files

    def __str__(self):
        return self.dot()

    def __len__(self):
        return self._len


class CandidateGroup:
    def __init__(self):
        self.id = str(uuid.uuid4())

    def __eq__(self, other):
        if isinstance(other, CandidateGroup):
            return self.id == other.id
        return False

    def __hash__(self):
        return hash(self.id)

    def merge(self, other):
        self.id = other.id
