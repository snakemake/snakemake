__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import asyncio
from builtins import ExceptionGroup
from collections import defaultdict
import os
import base64
from pathlib import Path
import tempfile
import json
import shutil
import functools

from itertools import chain, filterfalse
from operator import attrgetter
from typing import Iterable, List, Optional, Union
from collections.abc import AsyncGenerator
from abc import abstractmethod
from snakemake.settings.types import DeploymentMethod

from snakemake.template_rendering import check_template_output
from snakemake_interface_common.utils import lazy_property
from snakemake_interface_executor_plugins.jobs import (
    JobExecutorInterface,
    GroupJobExecutorInterface,
    SingleJobExecutorInterface,
)
from snakemake_interface_executor_plugins.settings import ExecMode
from snakemake_interface_logger_plugins.common import LogEvent

from snakemake.io import (
    _IOFile,
    IOFile,
    is_callable,
    Wildcards,
    Resources,
    is_flagged,
    get_flag_value,
    wait_for_files,
)
from snakemake.settings.types import SharedFSUsage
from snakemake.resources import GroupResources
from snakemake.target_jobs import TargetSpec
from snakemake.sourcecache import LocalSourceFile, SourceFile, infer_source_file
from snakemake.utils import format
from snakemake.exceptions import (
    InputOpenException,
    RuleException,
    ProtectedOutputException,
    WorkflowError,
)

from snakemake.logging import logger
from snakemake.common import (
    get_function_params,
    get_uuid,
    IO_PROP_LIMIT,
)
from snakemake.io.fmt import fmt_iofile
from snakemake.common.tbdstring import TBDString
from snakemake_interface_report_plugins.interfaces import JobReportInterface


def format_files(io, is_input: bool):
    return [fmt_iofile(f, as_input=is_input, as_output=not is_input) for f in io]


def jobfiles(jobs, type):
    return chain(*map(attrgetter(type), jobs))


@functools.lru_cache
def get_script_mtime(path: str) -> float:
    return os.lstat(path).st_mtime


class AbstractJob(JobExecutorInterface):
    @abstractmethod
    def reset_params_and_resources(self): ...

    @abstractmethod
    def get_target_spec(self): ...

    @abstractmethod
    def products(self, include_logfiles=True): ...

    def has_products(self, include_logfiles=True):
        for _ in self.products(include_logfiles=include_logfiles):
            return True
        return False

    def _get_scheduler_resources(self):
        if self._scheduler_resources is None:
            if self.dag.workflow.local_exec or self.is_local:
                res_dict = {
                    k: v
                    for k, v in self.resources.items()
                    if not isinstance(self.resources[k], TBDString)
                }
            else:
                res_dict = {
                    k: self.resources[k]
                    for k in (
                        set(self.resources.keys())
                        - self.dag.workflow.resource_scopes.locals
                    )
                    if not isinstance(self.resources[k], TBDString)
                }
            res_dict["_job_count"] = 1
            self._scheduler_resources = Resources(fromdict=res_dict)
        return self._scheduler_resources


class JobFactory:
    def __init__(self):
        self.cache = dict()

    async def new(
        self,
        rule,
        dag,
        wildcards_dict=None,
        format_wildcards=None,
        targetfile=None,
        update=False,
        groupid=None,
    ):
        key = (rule.name, *sorted(wildcards_dict.items()))

        async def new():
            if update:
                new_job = lambda: Job(
                    rule, dag, wildcards_dict, format_wildcards, targetfile, groupid
                )
            else:
                new_job = lambda: Job(
                    rule, dag, wildcards_dict, format_wildcards, targetfile
                )
            obj = None
            missing_iofiles = set()
            while obj is None:
                try:
                    obj = new_job()
                except InputOpenException as e:
                    if e.iofile in missing_iofiles:
                        raise WorkflowError(
                            f"Failed to open input file even after download from remote storage: {e.iofile}. "
                            "Has it been deleted by another process or is there a filesystem issue?",
                        )
                    if e.iofile.is_storage:
                        # retrieve the missing file from storage
                        await e.iofile.retrieve_from_storage()
                        missing_iofiles.add(e.iofile)
                    else:
                        raise WorkflowError(
                            f"Failed to open input file: {e.iofile}. Has it been deleted by another process?",
                            rule=e.rule,
                        )
            return obj

        if update:
            # cache entry has to be replaced because job shall be constructed from scratch
            obj = await new()
            self.cache[key] = obj
        else:
            try:
                # try to get job from cache
                obj = self.cache[key]
            except KeyError:
                obj = await new()
                self.cache[key] = obj
        return obj


class Job(AbstractJob, SingleJobExecutorInterface, JobReportInterface):
    obj_cache = dict()

    __slots__ = [
        "rule",
        "dag",
        "wildcards_dict",
        "_wildcards",
        "_format_wildcards",
        "_input",
        "dependencies",
        "_output",
        "_params",
        "_log",
        "_benchmark",
        "_resources",
        "_conda_env_file",
        "_conda_env",
        "_container_img_url",
        "_shadow_dir",
        "_inputsize",
        "temp_output",
        "protected_output",
        "touch_output",
        "_hash",
        "_attempt",
        "_group",
        "targetfile",
        "incomplete_input_expand",
        "_params_and_resources_resetted",
        "_queue_input",
        "_aux_resources",
        "_non_derived_params",
    ]

    def __init__(
        self,
        rule,
        dag,
        wildcards_dict=None,
        format_wildcards=None,
        targetfile=None,
        groupid=None,
    ):
        self.rule = rule
        self.dag = dag

        # the targetfile that led to the job
        # it is important to record this, since we need it to submit the
        # job on a cluster. In contrast, an arbitrary targetfile could
        # lead to a different composition of wildcard values (in case of
        # ambiguity in matching).
        self.targetfile = targetfile
        self.wildcards_dict = wildcards_dict
        self.wildcards = Wildcards(fromdict=self.wildcards_dict)
        self._format_wildcards = (
            self.wildcards
            if format_wildcards is None
            else Wildcards(fromdict=format_wildcards)
        )

        (
            self.input,
            input_mapping,
            self.dependencies,
            self.incomplete_input_expand,
        ) = self.rule.expand_input(self.wildcards_dict, groupid=groupid)

        self.output, output_mapping = self.rule.expand_output(self.wildcards_dict)
        # other properties are lazy to be able to use additional parameters and check already existing files
        self._params = None
        self._log = None
        self._benchmark = None
        self._resources = None
        self._conda_env_spec = None
        self._container_img_url = None
        self._scheduler_resources = None
        self._conda_env = None
        self._group = None

        # pipe_group will only be set if the job generates or consumes a pipe
        self.pipe_group = None

        self.shadow_dir = None
        self._inputsize = None
        self._is_updated = False
        self._params_and_resources_resetted = False
        self._non_derived_params = None

        self._attempt = self.dag.workflow.attempt
        self._aux_resources = dict()

        # TODO get rid of these
        self.temp_output, self.protected_output = set(), set()
        self.touch_output = set()
        self.pipe_or_service_output = set()
        self._queue_input = defaultdict(list)
        for f in self.output:
            f_ = output_mapping[f]
            if f_ in self.rule.temp_output:
                self.temp_output.add(f)
            if f_ in self.rule.protected_output:
                self.protected_output.add(f)
            if f_ in self.rule.touch_output:
                self.touch_output.add(f)
            if is_flagged(f_, "pipe") or is_flagged(f_, "service"):
                self.pipe_or_service_output.add(f)
        for f in self.input:
            f_ = input_mapping[f]
            queue_info = get_flag_value(f_, "from_queue")
            if queue_info:
                self._queue_input[queue_info].append(f)

    def add_aux_resource(self, name, value):
        if name in self._aux_resources:
            raise ValueError(
                f"Resource {name} already exists in aux_resources of job {self}."
            )
        self._aux_resources[name] = value

    @property
    def is_updated(self):
        return self._is_updated

    @is_updated.setter
    def is_updated(self, value):
        self._is_updated = value

    @property
    def shadow_dir(self):
        return self._shadow_dir

    @shadow_dir.setter
    def shadow_dir(self, value):
        self._shadow_dir = value

    @property
    def wildcards(self):
        return self._wildcards

    @wildcards.setter
    def wildcards(self, value):
        self._wildcards = value

    @property
    def input(self):
        return self._input

    @input.setter
    def input(self, value):
        self._input = value

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, value):
        self._output = value

    def logfile_suggestion(self, prefix: str) -> str:
        """Return a suggestion for the log file name given a prefix."""
        return (
            "/".join(
                [prefix, self.rule.name]
                + [
                    f"{w}_{v}"
                    for w, v in sorted(
                        self.wildcards_dict.items(), key=lambda item: item[0]
                    )
                ]
            )
            + ".log"
        )

    def has_queue_input(self) -> bool:
        return bool(self._queue_input)

    def has_unfinished_queue_input(self) -> bool:
        return any(not queue_info.finished for queue_info in self._queue_input)

    def queue_input_last_checked(self) -> float:
        return max(queue_info.last_checked or 0.0 for queue_info in self._queue_input)

    async def updated(self):
        group = self.dag.get_job_group(self)
        groupid = None
        if group is None:
            if self.dag.workflow.local_exec or self.is_local:
                groupid = self.dag.workflow.group_settings.local_groupid
        else:
            groupid = group.jobid

        job = await self.dag.job_factory.new(
            self.rule,
            self.dag,
            wildcards_dict=self.wildcards_dict,
            targetfile=self.targetfile,
            update=True,
            groupid=groupid,
        )
        job.is_updated = True
        return job

    def is_valid(self):
        """Check if job is valid"""
        # these properties have to work in dry-run as well. Hence we check them here:
        self.rule.expand_benchmark(self.wildcards_dict)
        self.rule.expand_log(self.wildcards_dict)

    async def outputs_older_than_script_or_notebook(self):
        """return output that's older than script, i.e. script has changed"""
        path = self.rule.script or self.rule.notebook
        if not path:
            return

        path: SourceFile = self._path_to_source_file(path)

        if isinstance(path, LocalSourceFile):
            path = path.get_path_or_uri()
            if os.path.exists(path):
                script_mtime = get_script_mtime(path)
                for f in self.output:
                    if await f.exists() and not await f.is_newer(script_mtime):
                        yield f
        # TODO also handle remote file case here.

    def get_target_spec(self):
        return [TargetSpec(self.rule.name, self.wildcards_dict)]

    @property
    def threads(self):
        return self.resources._cores

    @property
    def params(self):
        if self._params is None:
            self._expand_params()
        return self._params

    @property
    def non_derived_params(self):
        if self._non_derived_params is None:
            self._expand_params()
        return self._non_derived_params

    def _expand_params(self):
        self._params, self._non_derived_params = self.rule.expand_params(
            self.wildcards_dict, self.input, self.output, self
        )

    @property
    def log(self):
        if self._log is None:
            self._log = self.rule.expand_log(self.wildcards_dict)
        return self._log

    @property
    def benchmark(self):
        if self._benchmark is None:
            self._benchmark = self.rule.expand_benchmark(self.wildcards_dict)
        return self._benchmark

    @property
    def benchmark_repeats(self):
        if self.benchmark is not None:
            return get_flag_value(self.benchmark, "repeat") or 1

    @property
    def group(self):
        if self._group is None:
            self._group = self.rule.expand_group(self.wildcards_dict)
        return self._group

    @group.setter
    def group(self, group):
        self._group = group

    @property
    def attempt(self):
        return self._attempt

    @attempt.setter
    def attempt(self, attempt):
        # reset resources
        self._resources = None
        self._attempt = attempt

    def _get_resources_to_skip(self):
        """Return a set of resource names that are callable and depend on input files."""
        return {
            name
            for name, val in self.rule.resources.items()
            if is_callable(val) and "input" in get_function_params(val)
        }

    @property
    def resources(self):
        if self._resources is None:
            if self.dag.workflow.local_exec or self.is_local:
                skip_evaluation = set()
            else:
                # tmpdir should be evaluated in the context of the actual execution
                skip_evaluation = {"tmpdir"}
            if not self._params_and_resources_resetted:
                # initial evaluation, input files of job are probably not yet present.
                # Therefore skip all functions that depend on input files.
                skip_evaluation.update(self._get_resources_to_skip())
            self._resources = self.rule.expand_resources(
                self.wildcards_dict,
                self.input,
                self.attempt,
                skip_evaluation=skip_evaluation,
            )
        if self._aux_resources and any(
            name not in self._resources.keys() for name in self._aux_resources.keys()
        ):
            self._resources.update(self._aux_resources)
        return self._resources

    @property
    def scheduler_resources(self):
        return self._get_scheduler_resources()

    def reset_params_and_resources(self):
        if not self._params_and_resources_resetted:
            self._resources = None
            self._params = None
            self._params_and_resources_resetted = True

    @property
    def conda_env_spec(self):
        if self._conda_env_spec is None:
            self._conda_env_spec = self.rule.expand_conda_env(
                self.wildcards_dict, self.params, self.input
            )
        return self._conda_env_spec

    @property
    def conda_env(self):
        if self.conda_env_spec:
            if self._conda_env is None:
                self._conda_env = self.dag.conda_envs.get(
                    (self.conda_env_spec, self.container_img_url)
                )
            return self._conda_env
        return None

    def archive_conda_env(self):
        """Archive a conda environment into a custom local channel."""
        if self.conda_env_spec:
            if self.conda_env.is_externally_managed:
                raise WorkflowError(
                    "Workflow archives cannot be created for workflows using externally managed conda environments."
                    "Please use paths to YAML files for all your conda directives.",
                    rule=self.rule,
                )
            return self.conda_env.create_archive()
        return None

    @property
    def needs_singularity(self):
        return self.container_img is not None

    @property
    def container_img_url(self):
        if self._container_img_url is None:
            self._container_img_url = self.rule.expand_container_img(
                self.wildcards_dict
            )

        return self._container_img_url

    @property
    def is_containerized(self):
        return self.rule.is_containerized

    @property
    def container_img(self):
        if (
            DeploymentMethod.APPTAINER
            in self.dag.workflow.deployment_settings.deployment_method
            and self.container_img_url
        ):
            return self.dag.container_imgs[self.container_img_url]
        return None

    @property
    def env_modules(self):
        return self.rule.env_modules

    @property
    def container_img_path(self):
        return self.container_img.path if self.container_img else None

    @property
    def is_shadow(self):
        return self.rule.shadow_depth is not None

    @property
    def priority(self):
        return self.dag.priority(self)

    @property
    def b64id(self):
        return base64.b64encode(
            (self.rule.name + "".join(self.output)).encode("utf-8")
        ).decode("utf-8")

    async def inputsize(self):
        """
        Return the size of the input files.
        Input files need to be present.
        """
        if self._inputsize is None:
            self._inputsize = sum([await f.size() for f in self.input])
        return self._inputsize

    @property
    def message(self):
        """Return the message for this job."""
        try:
            return (
                self.format_wildcards(self.rule.message) if self.rule.message else None
            )
        except AttributeError as ex:
            raise RuleException(str(ex), rule=self.rule)
        except KeyError as ex:
            raise RuleException(
                "Unknown variable in message of shell command: {}".format(str(ex)),
                rule=self.rule,
            )

    @property
    def shellcmd(self):
        """Return the shell command."""
        try:
            return (
                self.format_wildcards(self.rule.shellcmd)
                if self.rule.shellcmd
                else None
            )
        except AttributeError as ex:
            raise RuleException(str(ex), rule=self.rule)
        except KeyError as ex:
            raise RuleException(
                "Unknown variable when printing shell command: {}".format(str(ex)),
                rule=self.rule,
            )

    def _path_to_source_file(self, path: Union[str, Path]) -> SourceFile:
        """Return the path to source file referred by rule, while expanding
        wildcards and params.
        """
        if isinstance(path, Path):
            path = str(path)

        if self.wildcards is not None or self.params is not None:
            # parse wildcards and params to get the correct script name
            path = format(path, wildcards=self.wildcards, params=self.params)
        path = infer_source_file(path)
        if isinstance(path, LocalSourceFile) and self.rule.basedir:
            # needed if rule is included from another subdirectory
            path = self.rule.basedir.join(path)

        return path

    @property
    def is_shell(self):
        return self.rule.is_shell

    @property
    def is_norun(self):
        return self.rule.norun

    @property
    def is_script(self):
        return self.rule.is_script

    @property
    def is_notebook(self):
        return self.rule.is_notebook

    @property
    def is_wrapper(self):
        return self.rule.is_wrapper

    @property
    def is_cwl(self):
        return self.rule.is_cwl

    @property
    def is_template_engine(self):
        return self.rule.is_template_engine

    @property
    def is_run(self):
        return not (
            self.is_shell
            or self.is_norun
            or self.is_script
            or self.is_notebook
            or self.is_wrapper
            or self.is_cwl
        )

    @property
    def is_pipe(self):
        return any(is_flagged(o, "pipe") for o in self.output)

    @property
    def is_service(self):
        return any(is_flagged(o, "service") for o in self.output)

    def shadowed_path(self, f):
        """Get the shadowed path of IOFile f."""
        if not self.shadow_dir:
            return f
        f_ = IOFile(os.path.join(self.shadow_dir, f), self.rule)
        # The shadowed path does not need the storage object, storage will be handled
        # after shadowing.
        f_.clone_flags(f, skip_storage_object=True)
        return f_

    @property
    async def output_mintime(self):
        """Return oldest output file."""

        async def get_mtime(f):
            if await f.exists():
                mtime = await f.mtime()
                return mtime.local_or_storage()
            else:
                return None

        async with asyncio.TaskGroup() as tg:
            tasks = [tg.create_task(get_mtime(f)) for f in self.output]
        mintimes = list(
            filter(lambda res: res is not None, (task.result() for task in tasks))
        )

        mintime = min(mintimes) if mintimes else None

        if self.benchmark and await self.benchmark.exists():
            mintime_benchmark = (await self.benchmark.mtime()).local_or_storage()
            if mintime is not None:
                return min(mintime, mintime_benchmark)
            else:
                return mintime_benchmark

        return mintime

    async def missing_output(self, requested):
        async def handle_file(f):
            # pipe or service output is always declared as missing
            # (even if it might be present on disk for some reason)
            if (
                is_flagged(f, "pipe")
                or is_flagged(f, "service")
                or not await f.exists()
            ):
                yield f

        for f in requested:
            async for f in handle_file(f):
                yield f

    @property
    def local_input(self):
        for f in self.input:
            if not f.is_storage:
                yield f

    @property
    def unique_input(self):
        seen = set()

        for element in filterfalse(seen.__contains__, self.input):
            seen.add(element)
            yield element

    @property
    def local_output(self):
        for f in self.output:
            if not f.is_storage:
                yield f

    @property
    def storage_input(self):
        for f in self.input:
            if f.is_storage:
                yield f

    @property
    def storage_output(self):
        for f in self.output:
            if f.is_storage:
                yield f

    async def existing_output(self) -> AsyncGenerator[_IOFile, None]:
        for f in self.output:
            if await f.exists():
                yield f

    async def check_protected_output(self):
        async with asyncio.TaskGroup() as tg:
            tasks = [tg.create_task(f.protected()) for f in self.output]
        protected = [f for f, task in zip(self.output, tasks) if task.result()]
        if protected:
            raise ProtectedOutputException(self, protected)

    async def remove_existing_output(self):
        """Clean up output before rules actually run"""
        for f, f_ in zip(self.output, self.rule.output):
            if is_flagged(f, "update"):
                # output files marked as to be updated are not removed
                continue
            try:
                # remove_non_empty_dir only applies to directories which aren't
                # flagged with directory().
                await f.remove(remove_non_empty_dir=False)
            except FileNotFoundError:
                # No file == no problem
                pass

        for f in self.log:
            await f.remove(remove_non_empty_dir=False)

    async def prepare(self):
        """
        Prepare execution of job.
        This includes creation of directories and deletion of previously
        created files.
        Creates a shadow directory for the job if specified.
        """

        await self.check_protected_output()

        unexpected_output = self.dag.reason(self).missing_output.intersection(
            [f async for f in self.existing_output()]
        )
        if unexpected_output:
            logger.warning(
                "Warning: the following output files of rule {} were not "
                "present when the DAG was created:\n{}".format(
                    self.rule, unexpected_output
                )
            )

        if not self.is_norun:
            await self.remove_existing_output()

            # Create tmpdir if necessary
            if self.resources.get("tmpdir"):
                os.makedirs(self.resources.tmpdir, exist_ok=True)

            for f, f_ in zip(self.output, self.rule.output):
                f.prepare()

            for f in self.log:
                f.prepare()
            if self.benchmark:
                self.benchmark.prepare()

        # wait for input files, respecting keep_storage_local
        try:
            await wait_for_files(
                self.input,
                wait_for_local=self.dag.workflow.keep_storage_local_at_runtime,
                latency_wait=self.dag.workflow.execution_settings.latency_wait,
                consider_local={
                    f for f in self.input if self.is_pipe_or_service_input(f)
                },
            )
        except IOError as ex:
            raise WorkflowError(
                ex,
                rule=self.rule,
            )

        if not self.is_shadow or self.is_norun:
            return

        # Create shadow directory structure
        self.shadow_dir = tempfile.mkdtemp(
            dir=self.rule.workflow.persistence.shadow_path
        )
        cwd = os.getcwd()

        # "minimal" creates symlinks only to the input files in the shadow directory
        # "copy-minimal" creates copies instead
        if (
            self.rule.shadow_depth == "minimal"
            or self.rule.shadow_depth == "copy-minimal"
        ):
            # Re-create the directory structure in the shadow directory
            for f, d in set(
                [
                    (item, os.path.dirname(item))
                    for sublist in [self.input, self.output, self.log]
                    if sublist is not None
                    for item in sublist
                ]
            ):
                if d and not os.path.isabs(d):
                    rel_path = os.path.relpath(d)
                    # Only create subdirectories
                    if rel_path.split(os.path.sep)[0] != "..":
                        os.makedirs(
                            os.path.join(self.shadow_dir, rel_path), exist_ok=True
                        )
                    else:
                        raise RuleException(
                            "The following file name references a parent directory relative to your workdir.\n"
                            'This isn\'t supported for shadow: "{}". Consider using an absolute path instead.\n{}'.format(
                                f, self.rule.shadow_depth
                            ),
                            rule=self.rule,
                        )

            # Symlink or copy the input files
            if self.rule.shadow_depth == "copy-minimal":
                for rel_path in set(
                    [os.path.relpath(f) for f in self.input if not os.path.isabs(f)]
                ):
                    copy = os.path.join(self.shadow_dir, rel_path)
                    shutil.copy(rel_path, copy)
            else:
                for rel_path in set(
                    [os.path.relpath(f) for f in self.input if not os.path.isabs(f)]
                ):
                    link = os.path.join(self.shadow_dir, rel_path)
                    original = os.path.relpath(rel_path, os.path.dirname(link))
                    os.symlink(original, link)

        # Shallow simply symlink everything in the working directory.
        elif self.rule.shadow_depth == "shallow":
            for source in os.listdir(cwd):
                link = os.path.join(self.shadow_dir, source)
                os.symlink(os.path.abspath(source), link)
        elif self.rule.shadow_depth == "full":
            snakemake_dir = os.path.join(cwd, ".snakemake")
            for dirpath, dirnames, filenames in os.walk(cwd, followlinks=True):
                # a link should not point to a parent directory of itself, else can cause infinite recursion
                # Must exclude .snakemake and its children to avoid infinite
                # loop of symlinks.
                if os.path.commonprefix([snakemake_dir, dirpath]) == snakemake_dir:
                    continue
                for dirname in dirnames:
                    if dirname == ".snakemake":
                        continue
                    relative_source = os.path.relpath(os.path.join(dirpath, dirname))
                    shadow = os.path.join(self.shadow_dir, relative_source)
                    os.mkdir(shadow)

                for filename in filenames:
                    source = os.path.join(dirpath, filename)
                    relative_source = os.path.relpath(source)
                    link = os.path.join(self.shadow_dir, relative_source)
                    os.symlink(source, link)

    def is_pipe_or_service_input(self, path) -> bool:
        for dep, files in self.dag.dependencies[self].items():
            if path in files:
                return path in dep.pipe_or_service_output
        return False

    async def cleanup(self):
        """Cleanup output files."""
        to_remove = [f for f in self.output if await f.exists()]
        to_remove.extend(
            [
                f
                for f in self.storage_output
                if (
                    await f.exists_in_storage()
                    if (f.is_storage and f.should_not_be_retrieved_from_storage)
                    else await f.exists_local()
                )
            ]
        )
        if to_remove:
            formatted = ", ".join(map(fmt_iofile, to_remove))
            logger.info(
                f"Removing output files of failed job {self}"
                f" since they might be corrupted:\n{formatted}"
            )
            for f in to_remove:
                await f.remove()

    def format_wildcards(self, string, **variables):
        """Format a string with variables from the job."""
        _variables = dict()
        _variables.update(self.rule.workflow.globals)
        _variables.update(
            dict(
                input=self.input,
                output=self.output,
                params=self.params,
                wildcards=self._format_wildcards,
                threads=self.threads,
                resources=self.resources,
                log=self.log,
                jobid=self.jobid,
                name=self.name,
                rule=self.rule.name,
                rulename=self.rule.name,
                bench_iteration=None,
            )
        )
        _variables.update(variables)
        try:
            return format(string, **_variables)
        except Exception as ex:
            raise RuleException(
                f"{ex.__class__.__name__}: {ex}, when formatting the following:\n"
                + string,
                rule=self.rule,
            )

    def properties(self, omit_resources=["_cores", "_nodes"], **aux_properties):
        resources = {
            name: res
            for name, res in self.resources.items()
            if name not in omit_resources
        }
        params = {name: value for name, value in self.params.items()}
        properties = {
            "type": "single",
            "rule": self.rule.name,
            "local": self.is_local,
            "input": None if len(self.input) > IO_PROP_LIMIT else self.input,
            "output": None if len(self.output) > IO_PROP_LIMIT else self.output,
            "wildcards": self.wildcards_dict,
            "params": params,
            "log": self.log,
            "threads": self.threads,
            "resources": resources,
            "jobid": self.dag.jobid(self),
        }
        properties.update(aux_properties)

        try:
            return json.dumps(properties)
        except TypeError:
            del properties["params"]
            return json.dumps(properties)

    @property
    def is_local(self):
        no_shared_fs = (
            SharedFSUsage.INPUT_OUTPUT
            not in self.dag.workflow.storage_settings.shared_fs_usage
        )
        return self.dag.workflow.is_local(self.rule) or (
            no_shared_fs
            and any(is_flagged(f, "local") for f in chain(self.input, self.output))
        )

    def __repr__(self):
        return self.rule.name

    def __lt__(self, other):
        return self.rule.__lt__(other.rule)

    def __gt__(self, other):
        return self.rule.__gt__(other.rule)

    def is_group(self):
        return False

    def log_info(self, indent=False, printshellcmd=True):
        priority = self.priority

        benchmark = (
            fmt_iofile(self.benchmark, as_input=False, as_output=True)
            if self.benchmark is not None
            else None
        )
        logger.info(
            f" Rule: {self.rule.name}, Jobid: {self.dag.jobid(self)}",
            extra=dict(
                event=LogEvent.JOB_INFO,
                jobid=self.dag.jobid(self),
                rule_msg=self.message,
                rule_name=self.rule.name,
                # in dryrun, we don't want to display a decision whether local or not
                # since we don't know how the user wants to execute
                local=(
                    not self.dag.workflow.dryrun
                    and self.dag.workflow.is_local(self.rule)
                ),
                input=format_files(self.input, is_input=True),
                output=format_files(self.output, is_input=False),
                log=format_files(self.log, is_input=False),
                benchmark=benchmark,
                wildcards=self.wildcards_dict,
                reason=str(self.dag.reason(self)),
                resources=self.resources,
                priority=(
                    "highest"
                    if priority == JobExecutorInterface.HIGHEST_PRIORITY
                    else priority
                ),
                threads=self.threads,
                indent=indent,
                is_checkpoint=self.rule.is_checkpoint,
                printshellcmd=printshellcmd,
                is_handover=self.rule.is_handover,
                shellcmd=self.shellcmd,
            ),
        )
        logger.info(
            f"Shell command: {self.shellcmd}",
            extra=dict(event=LogEvent.SHELLCMD, shellcmd=self.shellcmd, indent=indent),
        )
        if self.rule.is_checkpoint:
            logger.info("DAG of jobs will be updated after completion.")
        if self.rule.is_handover:
            logger.info("Handing over execution to foreign system...")

    def get_log_error_info(
        self, msg=None, indent=False, aux_logs: Optional[list] = None, **kwargs
    ):
        aux_logs = aux_logs or []
        # Retrieve conda env path only when conda is enabled with sdm
        # Otherwise the class Conda will be created also when not explicitly
        # requested with sdm, resulting in an error when conda is not available
        # in the container.
        conda_env_adress = (
            self.conda_env.address
            if (
                DeploymentMethod.CONDA
                in self.dag.workflow.deployment_settings.deployment_method
                and self.conda_env
            )
            else None
        )

        return dict(
            rule_name=self.rule.name,
            rule_msg=msg,
            jobid=self.dag.jobid(self),
            input=format_files(self.input, is_input=True),
            output=format_files(self.output, is_input=False),
            log=format_files(self.log, is_input=False) + aux_logs,
            conda_env=conda_env_adress,
            container_img=self.container_img,
            aux=kwargs,
            indent=indent,
            shellcmd=self.shellcmd,
            event=LogEvent.JOB_ERROR,
        )

    def log_error(
        self, msg=None, indent=False, aux_logs: Optional[list] = None, **kwargs
    ):
        logger.error(
            f"Error in rule {self.rule.name}, jobid: {self.dag.jobid(self)}",
            extra=self.get_log_error_info(msg, indent, aux_logs, **kwargs),
        )

    def register(self, external_jobid: Optional[str] = None):
        self.dag.workflow.persistence.started(self, external_jobid)

    def get_wait_for_files(self):
        wait_for_files = []
        wait_for_files.extend(self.local_input)
        wait_for_files.extend(
            f for f in self.storage_input if not f.should_not_be_retrieved_from_storage
        )

        if self.shadow_dir:
            wait_for_files.append(self.shadow_dir)
        if (
            DeploymentMethod.CONDA
            in self.dag.workflow.deployment_settings.deployment_method
            and self.conda_env
            and not self.conda_env.is_externally_managed
            and not self.conda_env.is_containerized
        ):
            # Managed or containerized envs are not present on the host FS,
            # hence we don't need to wait for them.
            wait_for_files.append(self.conda_env.address)
        return wait_for_files

    @property
    def jobid(self):
        return self.dag.jobid(self)

    def uuid(self):
        return str(
            get_uuid(
                f"{self.rule.name}:{','.join(sorted(f'{w}:{v}' for w, v in self.wildcards_dict.items()))}"
            )
        )

    async def postprocess(
        self,
        store_in_storage: bool = True,
        handle_log: bool = True,
        handle_touch: bool = True,
        error: bool = False,
        ignore_missing_output: Union[List[_IOFile], bool] = False,
        check_output_mtime: bool = True,
    ):
        if self.dag.is_draft_notebook_job(self):
            # no output produced but have to delete incomplete marker
            self.dag.workflow.persistence.cleanup(self)
            return

        shared_input_output = (
            SharedFSUsage.INPUT_OUTPUT
            in self.dag.workflow.storage_settings.shared_fs_usage
        )
        wait_for_local = (
            SharedFSUsage.STORAGE_LOCAL_COPIES
            in self.dag.workflow.storage_settings.shared_fs_usage
        )
        try:
            if (
                self.dag.workflow.exec_mode == ExecMode.SUBPROCESS
                or shared_input_output
                or (self.dag.workflow.remote_exec and not shared_input_output)
                or self.is_local
            ):
                if not error and handle_touch:
                    self.dag.handle_touch(self)
                if handle_log:
                    await self.dag.handle_log(self)
                if not error:
                    await self.dag.check_and_touch_output(
                        self,
                        wait=self.dag.workflow.execution_settings.latency_wait,
                        ignore_missing_output=ignore_missing_output,
                        # storage not yet handled, just require the local files
                        wait_for_local=True,
                        check_output_mtime=check_output_mtime,
                    )

                keep_shadow_dir = (
                    error and self.dag.workflow.execution_settings.keep_incomplete
                )
                kept_directory = self.dag.unshadow_output(
                    self, only_log=error, keep_shadow_dir=keep_shadow_dir
                )
                if kept_directory:
                    logger.error(
                        f"Keeping shadow directory: {kept_directory}.\n"
                        + "Run snakemake with --cleanup-shadow to clean up shadow directories."
                    )

                if (
                    not error
                    and self.rule.is_template_engine
                    and not is_flagged(self.output[0], "temp")
                ):
                    # TODO also check if consumers are executed on the same node
                    check_template_output(self)

                await self.dag.handle_storage(
                    self, store_in_storage=store_in_storage, store_only_log=error
                )
                if not error:
                    self.dag.handle_protected(self)
            elif not shared_input_output and not wait_for_local and not error:
                expanded_output = list(self.output)
                if self.benchmark:
                    expanded_output.append(self.benchmark)
                await wait_for_files(
                    expanded_output,
                    wait_for_local=False,
                    latency_wait=self.dag.workflow.execution_settings.latency_wait,
                    ignore_pipe_or_service=True,
                )

        except Exception as e:
            # cleanup metadata in case of any exception above
            self.dag.workflow.persistence.cleanup(self)
            raise e

        try:
            await self.dag.workflow.persistence.finished(self)
        except IOError as e:
            raise WorkflowError(
                "Error recording metadata for finished job "
                "({}). Please ensure write permissions for the "
                "directory {}".format(e, self.dag.workflow.persistence.path)
            )

        if error and not self.dag.workflow.execution_settings.keep_incomplete:
            await self.cleanup()

    @property
    def name(self):
        return self.rule.name

    @property
    def priority(self):
        return self.dag.priority(self)

    def products(self, include_logfiles=True):
        products = list(self.output)
        if self.benchmark:
            products.append(self.benchmark)
        if include_logfiles:
            products.extend(self.log)
        return products

    @property
    def rules(self):
        return [self.rule.name]

    @property
    def restart_times(self):
        return self.rule.restart_times

    @property
    def is_checkpoint(self):
        return self.rule.is_checkpoint

    def __len__(self):
        return 1


class GroupJobFactory:
    def __init__(self):
        self.cache = dict()

    def new(self, id, jobs, resources):
        jobs = frozenset(jobs)
        key = (id, jobs)
        try:
            obj = self.cache[key]
        except KeyError:
            obj = GroupJob(id, jobs, resources)
            self.cache[key] = obj
        return obj


class GroupJob(AbstractJob, GroupJobExecutorInterface):
    obj_cache = dict()

    __slots__ = [
        "_groupid",
        "_jobs",
        "_resources",
        "_input",
        "_output",
        "_log",
        "_inputsize",
        "_all_products",
        "_attempt",
        "_toposorted",
        "_jobid",
    ]

    def __init__(self, id, jobs, global_resources):
        self.groupid = id
        self._jobs = jobs
        self.global_resources = global_resources
        self._toposorted = None
        self._resources = None
        self._scheduler_resources = None
        self._input = None
        self._output = None
        self._log = None
        self._inputsize = None
        self._all_products = None
        self._attempt = self.dag.workflow.execution_settings.attempt
        self._jobid = None

    @property
    def groupid(self):
        return self._groupid

    @groupid.setter
    def groupid(self, new_groupid):
        self._groupid = new_groupid

    @property
    def jobs(self):
        return self._jobs

    @jobs.setter
    def jobs(self, new_jobs):
        self._jobs = new_jobs

    @property
    def is_containerized(self):
        return any(job.is_containerized for job in self.jobs)

    @property
    def toposorted(self):
        return self._toposorted

    @toposorted.setter
    def toposorted(self, new_toposorted):
        self._toposorted = new_toposorted

    def logfile_suggestion(self, prefix: str) -> str:
        """Return a suggestion for the log file name given a prefix."""
        return f"{prefix}/groupjobs/group_{self.name}/job_{self.jobid}.log"

    @property
    def dag(self):
        return next(iter(self.jobs)).dag

    async def retrieve_storage_input(self):
        for job in self.jobs:
            await job.retrieve_storage_input()

    def merge(self, other):
        assert other.groupid == self.groupid
        self.jobs = self.jobs | other.jobs

    def finalize(self):
        if self.toposorted is None:
            self.toposorted = [
                *self.dag.toposorted(self.jobs, inherit_pipe_dependencies=True)
            ]

    def __iter__(self):
        if self.toposorted is None:
            yield from self.jobs
            return

        for layer in self.toposorted:
            yield from layer

    def __repr__(self):
        return f"JobGroup({self.groupid},{repr(self.jobs)})"

    def __contains__(self, job):
        return job in self.jobs

    def is_group(self):
        return True

    @property
    def all_products(self):
        if self._all_products is None:
            self._all_products = set(f for job in self.jobs for f in job.products())
        return self._all_products

    @property
    def is_checkpoint(self):
        return any(job.is_checkpoint for job in self.jobs)

    @property
    def is_updated(self):
        return any(job.is_updated for job in self.jobs)

    def log_info(self):
        logger.info(
            f"Group job {self.groupid} (jobs in lexicogr. order):",
            extra=dict(
                event=LogEvent.GROUP_INFO, group_id=self.groupid, jobs=self.jobs
            ),
        )
        for job in sorted(self.jobs, key=lambda j: j.rule.name):
            job.log_info(indent=True)

    def log_error(self, msg=None, aux_logs: Optional[list] = None, **kwargs):
        job_error_info = [
            job.get_log_error_info(indent=True, **kwargs) for job in self.jobs
        ]
        aux_logs = aux_logs or []
        logger.error(
            f"Error in group {self.groupid}",
            dict(
                event=LogEvent.GROUP_ERROR,
                groupid=self.groupid,
                aux_logs=aux_logs,
                job_error_info=job_error_info,
                **kwargs,
            ),
        )
        for job in sorted(self.jobs, key=lambda j: j.rule.name):
            job.log_error(indent=True)

    def register(self, external_jobid: Optional[str] = None):
        for job in self.jobs:
            job.register(external_jobid=external_jobid)

    async def remove_existing_output(self):
        for job in self.jobs:
            await job.remove_existing_output()

    def reset_params_and_resources(self):
        for job in self.jobs:
            job.reset_params_and_resources()

    def get_wait_for_files(self):
        local_input = [
            f
            for job in self.jobs
            for f in job.local_input
            if f not in self.all_products
        ]
        remote_input = [
            f
            for job in self.jobs
            for f in job.storage_input
            if f not in self.all_products
        ]

        wait_for_files = []
        wait_for_files.extend(local_input)
        wait_for_files.extend(
            f for f in remote_input if not f.should_not_be_retrieved_from_storage
        )

        for job in self.jobs:
            if job.shadow_dir:
                wait_for_files.append(job.shadow_dir)
            if (
                DeploymentMethod.CONDA
                in self.dag.workflow.deployment_settings.deployment_method
                and job.conda_env
                and not job.conda_env.is_externally_managed
            ):
                wait_for_files.append(job.conda_env.address)
        return wait_for_files

    @property
    def resources(self):
        if self._resources is None:
            try:
                self._resources = GroupResources.basic_layered(
                    toposorted_jobs=self.toposorted,
                    constraints=self.global_resources,
                    run_local=self.dag.workflow.local_exec,
                    additive_resources=["runtime"],
                    sortby=["runtime"],
                )
            except WorkflowError as err:
                raise WorkflowError(
                    f"Error grouping resources in group '{self.groupid}': {err.args[0]}"
                )
        return Resources(fromdict=self._resources)

    @property
    def scheduler_resources(self):
        return self._get_scheduler_resources()

    @property
    def input(self):
        if self._input is None:
            self._input = [
                f for job in self.jobs for f in job.input if f not in self.all_products
            ]
        return self._input

    @property
    def output(self):
        all_input = set(f for job in self.jobs for f in job.input)
        if self._output is None:
            self._output = [
                f for job in self.jobs for f in job.output if f not in all_input
            ]
        return self._output

    @property
    def log(self):
        if self._log is None:
            self._log = [f for job in self.jobs for f in job.log]
        return self._log

    def products(self, include_logfiles=True):
        all_input = set(f for job in self.jobs for f in job.input)
        return [
            f
            for job in self.jobs
            for f in job.products(include_logfiles=include_logfiles)
            if f not in all_input
        ]

    def properties(self, omit_resources=["_cores", "_nodes"], **aux_properties):
        resources = {
            name: res
            for name, res in self.resources.items()
            if name not in omit_resources
        }
        properties = {
            "type": "group",
            "groupid": self.groupid,
            "local": self.is_local,
            "input": None if len(self.input) > IO_PROP_LIMIT else self.input,
            "output": None if len(self.output) > IO_PROP_LIMIT else self.output,
            "threads": self.threads,
            "resources": resources,
            "jobid": self.jobid,
        }
        properties.update(aux_properties)

        return json.dumps(properties)

    @property
    def jobid(self):
        if not self._jobid:
            # The uuid of the last job is sufficient to uniquely identify the group job.
            # This is true because each job can only occur in one group job.
            # Additionally, this is the most stable id we can get, even if the group
            # changes by adding more upstream jobs, e.g. due to groupid usage in input
            # functions (see Dag.update_incomplete_input_expand_jobs())
            last_job = sorted(self.toposorted[-1])[-1]
            self._jobid = last_job.uuid()
        return self._jobid

    async def cleanup(self):
        for job in self.jobs:
            await job.cleanup()

    async def postprocess(self, error=False, **kwargs):
        def needed(job_, f):
            # Files that are targetfiles, non-temp files, or temp files that
            # are needed by other unfinished jobs outside of this group
            # are still needed. Other files can be missing.
            return (
                f in self.dag.targetfiles
                or not is_flagged(f, "temp")
                or any(
                    f in files
                    for j, files in self.dag.depending[job_].items()
                    if (
                        not self.dag.finished(j)
                        and self.dag.needrun(j)
                        and j not in self.jobs
                    )
                )
            )

        # ignore missing output if the output is not needed outside of this group
        if not kwargs.get("ignore_missing_output", False):
            kwargs["ignore_missing_output"] = [
                f
                for j in self.jobs
                for f in j.output
                if is_flagged(f, "nodelocal")
                or (is_flagged(f, "temp") and not needed(j, f))
            ]

        # Iterate over jobs in toposorted order (see self.__iter__) to
        # ensure that outputs are touched in correct order.
        for i, level in enumerate(self.toposorted):
            if i > 0:
                # Wait a bit to ensure that subsequent levels really have distinct
                # modification times.
                await asyncio.sleep(0.1)
            try:
                async with asyncio.TaskGroup() as tg:
                    for job in level:
                        # postprocessing involves touching output files (to ensure that
                        # modification times are always correct. This has to happen in
                        # topological order, such that they are not mixed up.
                        # We have to disable output mtime checks here
                        # (at least for now), because they interfere with the update of
                        # intermediate file mtimes that might happen in previous levels.
                        tg.create_task(
                            job.postprocess(
                                error=error, check_output_mtime=False, **kwargs
                            )
                        )
            except ExceptionGroup as e:
                raise WorkflowError(
                    f"Error postprocessing group job {self.jobid}.",
                    *e.exceptions,
                )
        # remove all pipe and service outputs since all jobs of this group are done and the
        # outputs are no longer needed
        for job in self.jobs:
            for f in job.output:
                if is_flagged(f, "pipe") or is_flagged(f, "service"):
                    await f.remove()

    @property
    def name(self):
        rules = sorted({job.rule.name for job in self.jobs})
        if len(rules) > 5:
            rules = rules[:5] + ["..."]
        return f"{self.groupid}_{'_'.join(rules)}"

    def check_protected_output(self):
        for job in self.jobs:
            job.check_protected_output()

    async def inputsize(self):
        if self._inputsize is None:
            self._inputsize = sum([await f.size() for f in self.input])
        return self._inputsize

    @property
    def priority(self):
        return max(self.dag.priority(job) for job in self.jobs)

    @property
    def is_local(self):
        return any(job.is_local for job in self.jobs)

    def merged_wildcards(self):
        jobs = iter(self.jobs)
        merged_wildcards = Wildcards(toclone=next(jobs).wildcards)
        for job in jobs:
            for name, value in job.wildcards.items():
                if name not in merged_wildcards.keys():
                    merged_wildcards.append(value)
                    merged_wildcards._add_name(name)
        return merged_wildcards

    def format_wildcards(self, string, **variables):
        """Format a string with variables from the job."""

        _variables = dict()
        _variables.update(self.dag.workflow.globals)
        _variables.update(
            dict(
                input=self.input,
                output=self.output,
                threads=self.threads,
                wildcards=self.merged_wildcards(),
                jobid=self.jobid,
                name=self.name,
                rule="GROUP",
                rulename="GROUP",
                resources=self.resources,
            )
        )
        _variables.update(variables)
        try:
            return format(string, **_variables)
        except NameError as ex:
            raise WorkflowError(f"NameError with group job {self.jobid}: {str(ex)}")
        except IndexError as ex:
            raise WorkflowError(f"IndexError with group job {self.jobid}: {str(ex)}")
        except Exception as ex:
            raise WorkflowError(
                f"Error when formatting {string} for group job {self.jobid}: {ex}"
            )

    @property
    def threads(self):
        return self.resources["_cores"]

    def get_target_spec(self):
        return [TargetSpec(job.rule.name, job.wildcards_dict) for job in self.jobs]

    @property
    def attempt(self):
        return self._attempt

    @attempt.setter
    def attempt(self, attempt):
        # reset resources
        self._resources = None
        for job in self.jobs:
            job.attempt = attempt
        self._attempt = attempt

    @property
    def needs_singularity(self):
        return any(job.needs_singularity for job in self.jobs)

    @property
    def rules(self):
        return set(job.rule.name for job in self.jobs)

    @property
    def expanded_output(self):
        """Yields the entire expanded output of all jobs"""
        for job in self.jobs:
            yield from job.output

    @property
    def restart_times(self):
        return max(job.restart_times for job in self.jobs)

    def __len__(self):
        return len(self.jobs)

    def __hash__(self):
        return hash(self.jobs)

    def __eq__(self, other):
        if not isinstance(other, AbstractJob):
            return False
        if other.is_group():
            return self.jobs == other.jobs
        else:
            return False


class Reason:
    __slots__ = [
        "_updated_input",
        "_updated_input_run",
        "_missing_output",
        "_incomplete_output",
        "input_changed",
        "code_changed",
        "params_changed",
        "software_stack_changed",
        "forced",
        "noio",
        "nooutput",
        "unfinished_queue_input",
        "derived",
        "pipe",
        "service",
        "target",
        "finished",
        "cleanup_metadata_instructions",
        "no_metadata",
        "outdated_metadata",
    ]

    def __init__(self):
        from snakemake.persistence import NO_PARAMS_CHANGE

        self.finished = False
        self._updated_input = None
        self._updated_input_run = None
        self._missing_output = None
        self._incomplete_output = None
        self.params_changed = NO_PARAMS_CHANGE
        self.code_changed = False
        self.software_stack_changed = False
        self.input_changed = False
        self.forced = False
        self.noio = False
        self.nooutput = False
        self.derived = True
        self.pipe = False
        self.service = False
        self.cleanup_metadata_instructions = None
        self.unfinished_queue_input = False
        self.no_metadata = False
        self.outdated_metadata = False

    def set_cleanup_metadata_instructions(self, job):
        self.cleanup_metadata_instructions = (
            "To ignore these changes, run snakemake "
            f"--cleanup-metadata {' '.join(job.output)}"
        )

    def is_provenance_triggered(self):
        """Return True if reason is triggered by provenance information."""
        return (
            self.params_changed
            or self.code_changed
            or self.software_stack_changed
            or self.input_changed
        )

    @lazy_property
    def updated_input(self):
        return set()

    @lazy_property
    def updated_input_run(self):
        return set()

    @lazy_property
    def missing_output(self):
        return set()

    @lazy_property
    def incomplete_output(self):
        return set()

    def mark_finished(self):
        "called if the job has been run"
        self.finished = True

    def get_names(self):
        if self.forced:
            yield "forced"
        if self.noio:
            yield "neither input nor output"
        if self.nooutput:
            yield "run or shell but no output"
        if self._missing_output:
            yield "output files have to be generated"
        if self._incomplete_output:
            yield "incomplete output files"
        if self._updated_input:
            yield "updated input files"
        if self._updated_input_run:
            yield "input files updated by another job"
        if self.pipe:
            yield "pipe output needed by consuming job"
        if self.service:
            yield "provides service for consuming job"
        if self.input_changed:
            yield "set of input files has changed since last execution"
        if self.unfinished_queue_input:
            yield "input is yielded by queue that has not been terminated yet"
        if self.code_changed:
            yield "code has changed since last execution"
        if self.params_changed:
            yield "params have changed since last execution"
        if self.software_stack_changed:
            yield "software environment definition has changed since last execution"

    def __str__(self):
        def concat_files(files, is_input: bool):
            return ", ".join(format_files(files, is_input=is_input))

        s = list()
        if self.forced:
            s.append("Forced execution")
        else:
            if self.noio:
                s.append(
                    "Rules with neither input nor output files are always executed."
                )
            elif self.nooutput:
                s.append(
                    "Rules with a run or shell declaration but no output "
                    "are always executed."
                )
            else:
                if self._missing_output:
                    s.append(
                        f"Missing output files: {concat_files(self.missing_output, is_input=False)}"
                    )
                if self._incomplete_output:
                    s.append(
                        f"Incomplete output files: {concat_files(self.incomplete_output, is_input=False)}"
                    )
                if self._updated_input:
                    updated_input = self.updated_input - self.updated_input_run
                    s.append(
                        f"Updated input files: {concat_files(updated_input, is_input=True)}"
                    )
                if self._updated_input_run:
                    s.append(
                        f"Input files updated by another job: {concat_files(self.updated_input_run, is_input=True)}"
                    )
                if self.pipe:
                    s.append(
                        "Output file is a pipe and has to be filled for consuming job."
                    )
                if self.service:
                    s.append(
                        "Job provides a service which has to be kept active until all consumers are finished."
                    )

                if self.input_changed:
                    s.append("Set of input files has changed since last execution")
                if self.code_changed:
                    s.append("Code has changed since last execution")
                if self.params_changed:
                    s.append(
                        f"Params have changed since last execution: {self.params_changed}"
                    )
                if self.software_stack_changed:
                    s.append(
                        "Software environment definition has changed since last execution"
                    )

        s = "; ".join(s)
        if self.finished:
            return f"Finished (was: {s})"
        return s

    def __bool__(self):
        v = bool(
            self.updated_input
            or self.missing_output
            or self.forced
            or self.updated_input_run
            or self.noio
            or self.nooutput
            or self.pipe
            or self.service
            or self.code_changed
            or self.params_changed
            or self.software_stack_changed
            or self.input_changed
            or self.unfinished_queue_input
        )
        return v and not self.finished


def jobs_to_rulenames(jobs: Iterable[Job]) -> List[str]:
    return sorted({job.rule.name for job in jobs})
