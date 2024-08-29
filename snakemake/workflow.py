__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import copy
import hashlib
import os
import re
import subprocess
import sys
import tarfile
import tempfile
from collections import OrderedDict, namedtuple
from collections.abc import Mapping
from dataclasses import dataclass, field
from functools import partial, wraps
from itertools import chain, filterfalse
from pathlib import Path
import tarfile
import tempfile
from typing import Dict, Iterable, List, Optional, Set, Union
from .common.workdir_handler import WorkdirHandler
from .settings.types import (
    ConfigSettings,
    DAGSettings,
    DeploymentMethod,
    DeploymentSettings,
    ExecutionSettings,
    GroupSettings,
    OutputSettings,
    RemoteExecutionSettings,
    RerunTrigger,
    ResourceSettings,
    SchedulingSettings,
    StorageSettings,
    WorkflowSettings,
    SharedFSUsage,
)

from snakemake_interface_executor_plugins.workflow import WorkflowExecutorInterface
from snakemake_interface_executor_plugins.cli import (
    SpawnedJobArgsFactoryExecutorInterface,
)
from snakemake_interface_common.utils import lazy_property
from snakemake_interface_executor_plugins.settings import ExecutorSettingsBase
from snakemake_interface_executor_plugins.registry.plugin import (
    Plugin as ExecutorPlugin,
)
from snakemake_interface_executor_plugins.settings import ExecMode
from snakemake_interface_common.plugin_registry.plugin import TaggedSettings
from snakemake_interface_report_plugins.settings import ReportSettingsBase
from snakemake_interface_report_plugins.registry.plugin import Plugin as ReportPlugin

from .logging import logger, format_resources
from .rules import Rule, Ruleorder, RuleProxy
from .exceptions import (
    CreateCondaEnvironmentException,
    RuleException,
    CreateRuleException,
    UnknownRuleException,
    NoRulesException,
    WorkflowError,
    update_lineno,
)
from .dag import DAG, ChangeType
from .scheduler import JobScheduler
from .parser import parse
import snakemake.io
from .io import (
    protected,
    temp,
    temporary,
    ancient,
    directory,
    expand,
    glob_wildcards,
    flag,
    touch,
    unpack,
    local,
    pipe,
    service,
    repeat,
    report,
    multiext,
    ensure,
    from_queue,
    IOFile,
    sourcecache_entry,
)

from .persistence import Persistence
from .utils import update_config
from .script import script
from .notebook import notebook
from .wrapper import wrapper
from .cwl import cwl
from .template_rendering import render_template
from snakemake_interface_common.utils import not_iterable

import snakemake.wrapper
from .common import (
    ON_WINDOWS,
    async_run,
    get_appdirs,
    is_local_file,
    Rules,
    Scatter,
    Gather,
    smart_join,
    NOTHING_TO_BE_DONE_MSG,
)
from .utils import simplify_path
from .checkpoints import Checkpoints
from .resources import ResourceScopes
from .caching.local import OutputFileCache as LocalOutputFileCache
from .caching.storage import OutputFileCache as StorageOutputFileCache
from .modules import ModuleInfo, WorkflowModifier, get_name_modifier_func
from .ruleinfo import InOutput, RuleInfo
from .sourcecache import (
    LocalSourceFile,
    SourceCache,
    SourceFile,
    infer_source_file,
)
from .deployment.conda import Conda
from . import api, sourcecache
import snakemake.ioutils
import snakemake.ioflags
from .jobs import jobs_to_rulenames

SourceArchiveInfo = namedtuple("SourceArchiveInfo", ("query", "checksum"))


def rule_decorate(func):
    @wraps(func)
    def _(self: "Workflow", *n, **kw):
        def _d(ruleinfo: RuleInfo):
            func(self, ruleinfo, *n, **kw)
            return ruleinfo

        return _d

    return _


@dataclass
class Workflow(WorkflowExecutorInterface):
    config_settings: ConfigSettings
    resource_settings: ResourceSettings
    workflow_settings: WorkflowSettings
    storage_settings: Optional[StorageSettings] = None
    dag_settings: Optional[DAGSettings] = None
    execution_settings: Optional[ExecutionSettings] = None
    deployment_settings: Optional[DeploymentSettings] = None
    scheduling_settings: Optional[SchedulingSettings] = None
    output_settings: Optional[OutputSettings] = None
    remote_execution_settings: Optional[RemoteExecutionSettings] = None
    group_settings: Optional[GroupSettings] = None
    executor_settings: ExecutorSettingsBase = None
    storage_provider_settings: Optional[Mapping[str, TaggedSettings]] = None
    check_envvars: bool = True
    cache_rules: Dict[str, str | bool] = field(default_factory=dict)
    overwrite_workdir: Optional[str | Path] = None
    _workdir_handler: Optional[WorkdirHandler] = field(init=False, default=None)
    injected_conda_envs: List = field(default_factory=list)

    def __post_init__(self):
        """
        Create the controller.
        """
        from snakemake.storage import StorageRegistry

        self.global_resources: dict = dict(self.resource_settings.resources)
        self.global_resources["_cores"] = self.resource_settings.cores
        self.global_resources["_nodes"] = self.resource_settings.nodes

        self._rules: dict[str, Rule] = OrderedDict()
        self.default_target = None
        self._workdir_init = os.path.abspath(os.curdir)
        self._ruleorder = Ruleorder()
        self._localrules = set()
        self._linemaps = dict()
        self.rule_count = 0
        self.included = []
        self.included_stack: list[SourceFile] = []
        self._persistence: Optional[Persistence] = None
        self._dag: Optional[DAG] = None
        self._onsuccess = lambda log: None
        self._onerror = lambda log: None
        self._onstart = lambda log: None
        self._rulecount = 0
        self._parent_groupids = dict()
        self.global_container_img = None
        self.global_is_containerized = False
        self.configfiles = list(self.config_settings.configfiles)
        self.report_text = None
        # environment variables to pass to jobs
        # These are defined via the "envvars:" syntax in the Snakefile itself
        self._envvars = set()
        self._scatter = dict(self.resource_settings.overwrite_scatter)
        self._resource_scopes = ResourceScopes.defaults()
        self._resource_scopes.update(self.resource_settings.overwrite_resource_scopes)
        self.modules: dict[str, ModuleInfo] = dict()
        self._snakemake_tmp_dir = tempfile.TemporaryDirectory(prefix="snakemake")

        self._sourcecache = SourceCache(self.source_cache_path)

        self._scheduler = None
        self._spawned_job_general_args = None
        self._executor_plugin = None
        self._storage_registry = StorageRegistry(self)
        self._source_archive = None

        _globals = globals()
        from snakemake.shell import shell

        _globals["shell"] = shell
        _globals["workflow"] = self
        _globals["checkpoints"] = Checkpoints()
        _globals["scatter"] = Scatter()
        _globals["gather"] = Gather()
        _globals["github"] = sourcecache.GithubFile
        _globals["gitlab"] = sourcecache.GitlabFile
        _globals["gitfile"] = sourcecache.LocalGitFile
        _globals["storage"] = self._storage_registry
        snakemake.ioutils.register_in_globals(_globals)
        snakemake.ioflags.register_in_globals(_globals)
        _globals["from_queue"] = from_queue

        self.vanilla_globals = dict(_globals)
        self.modifier_stack = [WorkflowModifier(self, globals=_globals)]
        self._output_file_cache = None
        self.cache_rules = dict()

        self.globals["config"] = copy.deepcopy(self.config_settings.overwrite_config)

    @property
    def parent_groupids(self):
        return self._parent_groupids

    def tear_down(self):
        for conda_env in self.injected_conda_envs:
            conda_env.deactivate()
        if self._workdir_handler is not None:
            self._workdir_handler.change_back()
        self._snakemake_tmp_dir.cleanup()

    @property
    def is_main_process(self):
        return self.exec_mode == ExecMode.DEFAULT

    @property
    def snakemake_tmp_dir(self) -> Path:
        return Path(self._snakemake_tmp_dir.name)

    def register_resource(self, name: str, value: Union[int, str]):
        self.global_resources[name] = value
        if self.scheduler is not None:
            # update the scheduler if it is already active
            self.scheduler.resources[name] = value

    @property
    def source_cache_path(self) -> Path:
        assert self.storage_settings is not None
        if SharedFSUsage.SOURCE_CACHE not in self.storage_settings.shared_fs_usage:
            return self.snakemake_tmp_dir / "source-cache"
        else:
            return Path(
                os.path.join(get_appdirs().user_cache_dir, "snakemake/source-cache")
            )

    @property
    def storage_registry(self):
        return self._storage_registry

    @property
    def source_archive(self):
        assert self._source_archive is not None, (
            "bug: source archive info accessed but source archive has not been "
            "uploaded to default storage provider before"
        )
        return self._source_archive

    def upload_sources(self):
        assert self.storage_settings is not None
        with tempfile.NamedTemporaryFile(suffix="snakemake-sources.tar.xz") as tf:
            self.write_source_archive(Path(tf.name))
            tf.flush()
            with open(tf.name, "rb") as f:
                checksum = hashlib.file_digest(f, "sha256").hexdigest()

            prefix = self.storage_settings.default_storage_prefix
            if prefix:
                prefix = f"{prefix}/"
            query = f"{prefix}snakemake-workflow-sources.{checksum}.tar.xz"

            self._source_archive = SourceArchiveInfo(query, checksum)

            obj = self.storage_registry.default_storage_provider.object(query)
            obj.set_local_path(Path(tf.name))
            logger.info("Uploading source archive to storage provider...")
            async_run(obj.managed_store())

    def write_source_archive(self, path: Path):
        def get_files():
            for f in self.dag.get_sources():
                if f.startswith(".."):
                    logger.warning(
                        "Ignoring source file {}. Only files relative "
                        "to the working directory are allowed.".format(f)
                    )
                    continue

                # The kubernetes API can't create secret files larger than 1MB.
                source_file_size = os.path.getsize(f)
                max_file_size = 10000000
                if source_file_size > max_file_size:
                    logger.warning(
                        "Skipping the source file for upload {f}. Its size "
                        "{source_file_size} exceeds "
                        "the maximum file size (10MB). Consider to provide the file as "
                        "input file instead.".format(
                            f=f, source_file_size=source_file_size
                        )
                    )
                    continue
                yield f

        assert path.suffixes == [".tar", ".xz"]
        with tarfile.open(path, "w:xz") as archive:
            for f in get_files():
                archive.add(f)

    @property
    def enable_cache(self):
        return (
            self.workflow_settings is not None
            and self.workflow_settings.cache is not None
        )

    def check_cache_rules(self):
        for rule in self.rules:
            cache_mode = self.cache_rules.get(rule.name)
            if cache_mode:
                if len(rule.output) > 1:
                    if not all(out.is_multiext for out in rule.output):
                        raise WorkflowError(
                            "Rule is marked for between workflow caching but has multiple output files. "
                            "This is only allowed if multiext() is used to declare them (see docs on between "
                            "workflow caching).",
                            rule=rule,
                        )
                if not self.enable_cache:
                    logger.warning(
                        f"Workflow defines that rule {rule.name} is eligible for caching between workflows "
                        "(use the --cache argument to enable this)."
                    )
                if rule.benchmark:
                    raise WorkflowError(
                        "Rules with a benchmark directive may not be marked as eligible "
                        "for between-workflow caching at the same time. The reason is that "
                        "when the result is taken from cache, there is no way to fill the benchmark file with "
                        "any reasonable values. Either remove the benchmark directive or disable "
                        "between-workflow caching for this rule.",
                        rule=rule,
                    )

    @property
    def attempt(self):
        if self.execution_settings is None:
            # if not executing, we can safely set this to 1
            return 1
        return self.execution_settings.attempt

    @property
    def executor_plugin(self):
        return self._executor_plugin

    @property
    def dryrun(self):
        if self.executor_plugin is None:
            return False
        else:
            return self.executor_plugin.common_settings.dryrun_exec

    @property
    def touch(self):
        import snakemake.executors.touch

        return issubclass(
            self.executor_plugin.executor, snakemake.executors.touch.Executor
        )

    @property
    def use_threads(self):
        assert self.execution_settings is not None
        return (
            (
                self.execution_settings is not None
                and self.execution_settings.use_threads
            )
            or (os.name not in ["posix", "nt"])
            or not self.local_exec
        )

    @property
    def local_exec(self):
        if self.executor_plugin is not None:
            return self.executor_plugin.common_settings.local_exec
        else:
            return True

    @property
    def non_local_exec(self):
        return not self.local_exec

    @property
    def remote_exec(self):
        return self.exec_mode == ExecMode.REMOTE

    @property
    def exec_mode(self):
        return self.workflow_settings.exec_mode

    @lazy_property
    def spawned_job_args_factory(self) -> SpawnedJobArgsFactoryExecutorInterface:
        from snakemake.spawn_jobs import SpawnedJobArgsFactory

        return SpawnedJobArgsFactory(self)

    @property
    def basedir(self):
        return os.path.dirname(self.main_snakefile)

    @property
    def scheduler(self):
        return self._scheduler

    @scheduler.setter
    def scheduler(self, scheduler):
        self._scheduler = scheduler

    @property
    def envvars(self):
        return self._envvars

    @property
    def sourcecache(self):
        return self._sourcecache

    @property
    def workdir_init(self):
        return self._workdir_init

    @property
    def linemaps(self):
        return self._linemaps

    @property
    def persistence(self):
        return self._persistence

    @property
    def dag(self):
        return self._dag

    @property
    def main_snakefile(self) -> str:
        return self.included[0].get_path_or_uri()

    @property
    def output_file_cache(self):
        return self._output_file_cache

    @property
    def resource_scopes(self):
        return self._resource_scopes

    @property
    def overwrite_configfiles(self):
        return self.config_settings.configfiles

    @property
    def rerun_triggers(self) -> Set[RerunTrigger]:
        assert self.dag_settings is not None
        return self.dag_settings.rerun_triggers  # type: ignore[return-value]

    @property
    def conda_base_path(self):
        assert self.deployment_settings is not None
        if self.deployment_settings.conda_base_path:
            return self.deployment_settings.conda_base_path
        if DeploymentMethod.CONDA in self.deployment_settings.deployment_method:
            try:
                return Conda().prefix_path
            except CreateCondaEnvironmentException:
                # Return no preset conda base path now and report error later in jobs.
                return None
        else:
            return None

    @property
    def modifier(self):
        return self.modifier_stack[-1]

    @property
    def wildcard_constraints(self):
        return self.modifier.wildcard_constraints

    @property
    def globals(self):
        return self.modifier.globals

    def lint(self, json=False):
        from snakemake.linting.rules import RuleLinter
        from snakemake.linting.snakefiles import SnakefileLinter

        json_snakefile_lints, snakefile_linted = SnakefileLinter(
            self, self.included
        ).lint(json=json)
        json_rule_lints, rules_linted = RuleLinter(self, self.rules).lint(json=json)

        linted = snakefile_linted or rules_linted

        if json:
            import json

            print(
                json.dumps(
                    {"snakefiles": json_snakefile_lints, "rules": json_rule_lints},
                    indent=2,
                )
            )
        else:
            if not linted:
                logger.info("Congratulations, your workflow is in a good condition!")
        return linted

    def get_cache_mode(self, rule: Rule):
        if self.workflow_settings.cache is None:
            return None
        else:
            return self.cache_rules.get(rule.name)

    @property
    def rules(self):
        return self._rules.values()

    @property
    def cores(self):
        if self._cores is None:
            raise WorkflowError(
                "Workflow requires a total number of cores to be defined (e.g. because a "
                "rule defines its number of threads as a fraction of a total number of cores). "
                "Please set it with --cores N with N being the desired number of cores. "
                "Consider to use this in combination with --max-threads to avoid "
                "jobs with too many threads for your setup. Also make sure to perform "
                "a dryrun first."
            )
        return self._cores

    @property
    def _cores(self):
        return self.global_resources["_cores"]

    @property
    def nodes(self):
        return self.global_resources["_nodes"]

    @property
    def concrete_files(self):
        return (
            file
            for rule in self.rules
            for file in chain(rule.input, rule.output)
            if not callable(file) and not file.contains_wildcard()
        )

    def check(self):
        for clause in self._ruleorder:
            for rulename in clause:
                if not self.is_rule(rulename):
                    raise UnknownRuleException(
                        rulename, prefix="Error in ruleorder definition."
                    )
        self.check_cache_rules()
        self.check_localrules()

    def add_rule(self, rule: Rule, checkpoint=False, allow_overwrite=False):
        """
        Add a rule.
        """
        is_overwrite = self.is_rule(rule.name)
        if not allow_overwrite and is_overwrite:
            raise CreateRuleException(
                f"The name {rule.name} is already used by another rule",
                lineno=rule.lineno,
                snakefile=rule.snakefile,
            )
        self._rules[rule.name] = rule
        self.modifier.rules.add(rule)
        if not is_overwrite:
            self.rule_count += 1
        if not self.default_target:
            self.default_target = rule.name
        return rule  # self.get_rule(name)

    def is_rule(self, name):
        """
        Return True if name is the name of a rule.

        Arguments
        name -- a name
        """
        return name in self._rules

    def get_rule(self, name):
        """
        Get rule by name.

        Arguments
        name -- the name of the rule
        """
        if not self._rules:
            raise NoRulesException()
        if name not in self._rules:
            raise UnknownRuleException(name)
        return self._rules[name]

    def list_rules(self, only_targets=False):
        rules = self.rules
        if only_targets:
            rules = filterfalse(Rule.has_wildcards, rules)
        for rule in sorted(rules, key=lambda r: r.name):
            docstring = f" ({rule.docstring})" if rule.docstring else ""
            print(rule.name + docstring)

    def list_resources(self):
        for resource in set(
            resource for rule in self.rules for resource in rule.resources
        ):
            if resource not in "_cores _nodes".split():
                logger.info(resource)

    def is_local(self, rule):
        return self.local_exec or (
            rule.group is None
            and (rule.name in self._localrules or rule.norun or rule.is_template_engine)
        )

    def check_localrules(self):
        undefined = self._localrules - set(rule.name for rule in self.rules)
        if undefined:
            logger.warning(
                "localrules directive specifies rules that are not "
                "present in the Snakefile:\n{}\n".format(
                    "\n".join(map("\t{}".format, undefined))
                )
            )

    def inputfile(self, path):
        """Mark file as being an input file of the workflow.

        This also means that eventual --default-remote-provider/prefix settings
        will be applied to this file. The file is returned as _IOFile object,
        such that it can e.g. be transparently opened with _IOFile.open().
        """
        if isinstance(path, Path):
            path = str(path)
        assert self.storage_settings is not None
        if self.storage_settings.default_storage_provider is not None:
            path = self.modifier.modify_path(path)
        return IOFile(path)

    def _prepare_dag(
        self,
        forceall: bool,
        ignore_incomplete: bool,
        lock_warn_only: bool,
        nolock: bool = False,
        shadow_prefix: str | Path | None = None,
    ):
        if self.workflow_settings.cache is not None:
            self.cache_rules.update(
                {rulename: "all" for rulename in self.workflow_settings.cache}
            )
            if (
                self.storage_settings is not None
                and self.storage_settings.default_storage_provider is not None
            ):
                self._output_file_cache = StorageOutputFileCache(
                    self.storage_registry.default_storage_provider
                )
            else:
                self._output_file_cache = LocalOutputFileCache()

        def rules(items):
            return map(self._rules.__getitem__, filter(self.is_rule, items))

        assert self.dag_settings is not None
        if self.dag_settings.target_files_omit_workdir_adjustment:

            def files(items):
                return map(
                    self.modifier.path_modifier.apply_default_storage,
                    filterfalse(self.is_rule, items),
                )

        else:

            def files(items):
                relpath = lambda f: (
                    f
                    if os.path.isabs(f) or f.startswith("root://")
                    else os.path.relpath(f)
                )
                return map(
                    self.modifier.path_modifier.apply_default_storage,
                    map(relpath, filterfalse(self.is_rule, items)),
                )

        self.iocache = snakemake.io.IOCache(self.dag_settings.max_inventory_wait_time)

        if not self.dag_settings.targets and not self.dag_settings.target_jobs:
            targets: Iterable = (
                [self.default_target] if self.default_target is not None else list()
            )
        else:
            targets = self.dag_settings.targets

        prioritytargets = (
            set()
            if self.scheduling_settings is None
            else self.scheduling_settings.prioritytargets
        )

        priorityrules = set(rules(prioritytargets))
        priorityfiles = set(files(prioritytargets))
        forcerules = set(rules(self.dag_settings.forcerun))
        forcefiles = set(files(self.dag_settings.forcerun))
        untilrules = set(rules(self.dag_settings.until))
        untilfiles = set(files(self.dag_settings.until))
        omitrules = set(rules(self.dag_settings.omit_from))
        omitfiles = set(files(self.dag_settings.omit_from))
        targetrules = set(
            chain(
                rules(targets),
                filterfalse(Rule.has_wildcards, priorityrules),
                filterfalse(Rule.has_wildcards, forcerules),
                filterfalse(Rule.has_wildcards, untilrules),
            )
        )
        targetfiles = set(chain(files(targets), priorityfiles, forcefiles, untilfiles))

        if ON_WINDOWS:
            targetfiles = set(tf.replace(os.sep, os.altsep) for tf in targetfiles)

        if self.dag_settings.forcetargets:
            forcefiles.update(targetfiles)
            forcerules.update(targetrules)

        rules = (
            [
                rule
                for rule in self.rules
                if rule.name in self.dag_settings.allowed_rules
            ]
            if self.dag_settings.allowed_rules
            else self.rules
        )

        self._dag = DAG(
            self,
            rules,
            targetfiles=targetfiles,
            targetrules=targetrules,
            # when cleaning up conda or containers, we should enforce all possible jobs
            # since their envs shall not be deleted
            forceall=forceall,
            forcefiles=forcefiles,
            forcerules=forcerules,
            priorityfiles=priorityfiles,
            priorityrules=priorityrules,
            untilfiles=untilfiles,
            untilrules=untilrules,
            omitfiles=omitfiles,
            omitrules=omitrules,
            ignore_incomplete=ignore_incomplete,
        )

        persistence_path = (
            self.snakemake_tmp_dir / "persistence"
            if (
                self.storage_settings is not None
                and SharedFSUsage.PERSISTENCE
                not in self.storage_settings.shared_fs_usage
            )
            else None
        )

        assert self.deployment_settings is not None
        self._persistence = Persistence(
            nolock=nolock,
            dag=self._dag,
            conda_prefix=self.deployment_settings.conda_prefix,
            singularity_prefix=self.deployment_settings.apptainer_prefix,
            shadow_prefix=shadow_prefix,
            warn_only=lock_warn_only,
            path=persistence_path,
        )

    def generate_unit_tests(self, path: Path):
        """Generate unit tests for the workflow.

        Arguments
        path -- Path to the directory where the unit tests shall be generated.
        """
        from snakemake import unit_tests

        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=False,
            lock_warn_only=False,
        )
        self._build_dag()

        deploy = []
        assert self.deployment_settings is not None
        if DeploymentMethod.CONDA in self.deployment_settings.deployment_method:
            deploy.append("conda")
        if DeploymentMethod.APPTAINER in self.deployment_settings.deployment_method:
            deploy.append("singularity")
        unit_tests.generate(
            self.dag, path, deploy, configfiles=self.overwrite_configfiles
        )

    def cleanup_metadata(self, paths: List[Path]):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=False,
        )
        failed = []
        for path in paths:
            success = self.persistence.cleanup_metadata(path)
            if not success:
                failed.append(str(path))
        if failed:
            raise WorkflowError(
                "Failed to clean up metadata for the following files because the metadata was not present.\n"
                "If this is expected, there is nothing to do.\nOtherwise, the reason might be file system latency "
                "or still running jobs.\nConsider running metadata cleanup again.\nFiles:\n"
                + "\n".join(failed)
            )

    def unlock(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=False,
        )
        self._build_dag()
        try:
            self.persistence.cleanup_locks()
            logger.info("Unlocked working directory.")
        except IOError as e:
            raise WorkflowError(
                f"Error: Unlocking the directory {os.getcwd()} failed. Maybe "
                "you don't have the permissions?",
                e,
            )

    def cleanup_shadow(self):
        self._prepare_dag(forceall=False, ignore_incomplete=False, lock_warn_only=False)
        self._build_dag()
        with self.persistence.lock():
            self.persistence.cleanup_shadow()

    def delete_output(self, only_temp: bool = False, dryrun: bool = False):
        self._prepare_dag(forceall=False, ignore_incomplete=False, lock_warn_only=True)
        self._build_dag()

        async_run(self.dag.clean(only_temp=only_temp, dryrun=dryrun))

    def list_untracked(self):
        self._prepare_dag(forceall=False, ignore_incomplete=False, lock_warn_only=True)
        self._build_dag()

        self.dag.list_untracked()

    def list_changes(self, change_type: ChangeType):
        self._prepare_dag(forceall=False, ignore_incomplete=False, lock_warn_only=True)
        self._build_dag()

        items = async_run(self.dag.get_outputs_with_changes(change_type))
        if items:
            print(*items, sep="\n")

    def archive(self, path: Path):
        """Archive the workflow.

        Arguments
        path -- Path to the archive file.
        """
        self._prepare_dag(forceall=False, ignore_incomplete=False, lock_warn_only=True)
        self._build_dag()

        self.dag.archive(path)

    def summary(self, detailed: bool = False):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=True,
        )
        self._build_dag()

        async def join_summary(detailed):
            return "\n".join(
                [line async for line in self.dag.summary(detailed=detailed)]
            )

        print(async_run(join_summary(detailed)))

    def printdag(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=True,
        )
        self._build_dag()
        print(self.dag)

    def printrulegraph(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=True,
        )
        self._build_dag()
        print(self.dag.rule_dot())

    def printfilegraph(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=True,
        )
        self._build_dag()
        print(self.dag.filegraph_dot())

    def printd3dag(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=True,
        )
        self._build_dag()

        self.dag.d3dag()

    def containerize(self):
        from snakemake.deployment.containerize import containerize

        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=False,
            lock_warn_only=False,
        )
        self._build_dag()
        with self.persistence.lock():
            containerize(self, self.dag)

    def export_cwl(self, path: Path):
        """Export the workflow as CWL document.

        Arguments
        path -- the path to the CWL document to be created.
        """
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=False,
        )
        self._build_dag()

        import json

        from snakemake.cwl import dag_to_cwl

        with open(path, "w") as cwl:
            json.dump(dag_to_cwl(self.dag), cwl, indent=4)

    def create_report(
        self, report_plugin: ReportPlugin, report_settings: ReportSettingsBase
    ):
        from snakemake.report import auto_report

        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=False,
            lock_warn_only=False,
        )
        self._build_dag()

        async_run(auto_report(self.dag, report_plugin, report_settings))

    def conda_list_envs(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=False,
            lock_warn_only=False,
        )
        self._build_dag()
        self.dag.create_conda_envs(
            dryrun=True,
            quiet=True,
        )
        print("environment", "container", "location", sep="\t")
        for env in set(job.conda_env for job in self.dag.jobs):
            if env and not env.is_named:
                print(
                    env.file.simplify_path(),
                    env.container_img_url or "",
                    simplify_path(env.address),
                    sep="\t",
                )
        return True

    def conda_create_envs(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=False,
        )
        self._build_dag()

        assert self.deployment_settings is not None
        if DeploymentMethod.APPTAINER in self.deployment_settings.deployment_method:
            self.dag.pull_container_imgs()
        self.dag.create_conda_envs()

    def conda_cleanup_envs(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,  # True?
            ignore_incomplete=True,
            lock_warn_only=False,
        )
        self._build_dag()
        self.persistence.conda_cleanup_envs()

    def container_cleanup_images(self):
        assert self.dag_settings is not None
        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=True,
            lock_warn_only=False,
        )
        self._build_dag()
        self.persistence.cleanup_containers()

    def _build_dag(self):
        logger.info("Building DAG of jobs...")
        async_run(self.dag.init())
        async_run(self.dag.update_checkpoint_dependencies())

    def execute(
        self,
        executor_plugin: ExecutorPlugin,
        executor_settings: ExecutorSettingsBase,
        updated_files: Optional[List[str]] = None,
    ):
        logger.host_info()

        from snakemake.shell import shell

        assert self.deployment_settings is not None
        assert self.execution_settings is not None
        assert self.storage_settings is not None
        assert self.dag_settings is not None
        assert self.remote_execution_settings is not None
        assert self.output_settings is not None
        shell.conda_block_conflicting_envvars = (
            not self.deployment_settings.conda_not_block_search_path_envvars
        )

        if self.remote_execution_settings.envvars:
            self.register_envvars(*self.remote_execution_settings.envvars)

        self._executor_plugin = executor_plugin
        self.executor_settings = executor_settings

        if self.execution_settings.wait_for_files:
            try:
                async_run(
                    snakemake.io.wait_for_files(
                        self.execution_settings.wait_for_files,
                        latency_wait=self.execution_settings.latency_wait,
                    )
                )
            except IOError as e:
                logger.error(str(e))
                return False

        self._prepare_dag(
            forceall=self.dag_settings.forceall,
            ignore_incomplete=self.execution_settings.ignore_incomplete,
            lock_warn_only=self.dryrun,
            nolock=not self.execution_settings.lock,
            shadow_prefix=self.execution_settings.shadow_prefix,
        )

        if self.exec_mode in [ExecMode.SUBPROCESS, ExecMode.REMOTE]:
            self.persistence.deactivate_cache()

        self._build_dag()

        with self.persistence.lock():
            async_run(self.dag.postprocess(update_needrun=False))
            if not self.dryrun:
                # deactivate IOCache such that from now on we always get updated
                # size, existence and mtime information
                # ATTENTION: this may never be removed without really good reason.
                # Otherwise weird things may happen.
                self.iocache.deactivate()
                # clear and deactivate persistence cache, from now on we want to see updates
                self.persistence.deactivate_cache()

            if self.remote_execution_settings.immediate_submit and any(
                self.dag.checkpoint_jobs
            ):
                raise WorkflowError(
                    "Immediate submit mode (--immediate-submit) may not be used for workflows "
                    "with checkpoint jobs, as the dependencies cannot be determined before "
                    "execution in such cases."
                )
            if self.touch:
                self.dag.check_touch_compatible()

            if updated_files is not None:
                updated_files.extend(
                    f for job in self.dag.needrun_jobs() for f in job.output
                )

            shared_deployment = (
                SharedFSUsage.SOFTWARE_DEPLOYMENT
                in self.storage_settings.shared_fs_usage
            )

            if shared_deployment or (self.remote_exec and not shared_deployment):
                if (
                    DeploymentMethod.APPTAINER
                    in self.deployment_settings.deployment_method
                ):
                    self.dag.pull_container_imgs()
                if DeploymentMethod.CONDA in self.deployment_settings.deployment_method:
                    self.dag.create_conda_envs()

            shared_storage_local_copies = (
                SharedFSUsage.STORAGE_LOCAL_COPIES
                in self.storage_settings.shared_fs_usage
            )
            logger.debug(f"shared_storage_local_copies: {shared_storage_local_copies}")
            logger.debug(f"remote_exec: {self.remote_exec}")
            dryrun_or_touch = self.dryrun or self.touch
            if not dryrun_or_touch and (
                (self.exec_mode == ExecMode.DEFAULT and shared_storage_local_copies)
                or (self.remote_exec and not shared_storage_local_copies)
            ):
                async_run(self.dag.retrieve_storage_inputs())

            if (
                SharedFSUsage.SOURCES not in self.storage_settings.shared_fs_usage
                and self.exec_mode == ExecMode.DEFAULT
                and self.remote_execution_settings.job_deploy_sources
                and not executor_plugin.common_settings.can_transfer_local_files
            ):
                # no shared FS, hence we have to upload the sources to the storage
                self.upload_sources()

            def log_missing_metadata_info():
                no_metadata_jobs = [
                    job
                    for job in self.dag.jobs
                    if self.dag.reason(job).no_metadata and not self.dag.needrun(job)
                ]
                if no_metadata_jobs:
                    logger.info(
                        f"{len(no_metadata_jobs)} jobs have no recorded "
                        "provenance/metadata so that they "
                        "cannot be triggered by that. \n"
                        "Rules with missing "
                        f"metadata: {' '.join(jobs_to_rulenames(no_metadata_jobs))}"
                    )

            self.scheduler = JobScheduler(self, executor_plugin)

            if not self.dryrun:
                if len(self.dag):
                    from snakemake.shell import shell

                    shell_exec = shell.get_executable()
                    if shell_exec is not None:
                        logger.info(f"Using shell: {shell_exec}")
                    if not self.local_exec:
                        logger.resources_info(f"Provided remote nodes: {self.nodes}")
                    else:
                        if self._cores is not None:
                            warning = (
                                ""
                                if self._cores > 1
                                else " (use --cores to define parallelism)"
                            )
                            logger.resources_info(
                                f"Provided cores: {self._cores}{warning}"
                            )
                            logger.resources_info(
                                "Rules claiming more threads will be scaled down."
                            )

                    provided_resources = format_resources(self.global_resources)
                    if provided_resources:
                        logger.resources_info(
                            f"Provided resources: {provided_resources}"
                        )

                    if self.local_exec and any(rule.group for rule in self.rules):
                        logger.info("Group jobs: inactive (local execution)")

                    if (
                        DeploymentMethod.CONDA
                        not in self.deployment_settings.deployment_method
                        and any(rule.conda_env for rule in self.rules)
                    ):
                        logger.info("Conda environments: ignored")

                    if (
                        DeploymentMethod.APPTAINER
                        not in self.deployment_settings.deployment_method
                        and any(rule.container_img for rule in self.rules)
                    ):
                        logger.info("Singularity containers: ignored")

                    if self.exec_mode == ExecMode.DEFAULT:
                        logger.run_info("\n".join(self.dag.stats()))
                else:
                    logger.info(NOTHING_TO_BE_DONE_MSG)
                    return
            else:
                # the dryrun case
                if len(self.dag):
                    logger.run_info("\n".join(self.dag.stats()))
                else:
                    logger.info(NOTHING_TO_BE_DONE_MSG)
                    log_missing_metadata_info()
                    return
                if self.output_settings.quiet:
                    # in case of dryrun and quiet, just print above info and exit
                    return

            if not self.dryrun and not self.execution_settings.no_hooks:
                self._onstart(logger.get_logfile())

            def log_provenance_info():
                provenance_triggered_jobs = [
                    job
                    for job in self.dag.needrun_jobs(exclude_finished=False)
                    if self.dag.reason(job).is_provenance_triggered()
                ]
                if provenance_triggered_jobs:
                    logger.info(
                        "Some jobs were triggered by provenance information, "
                        "see 'reason' section in the rule displays above.\n"
                        "If you prefer that only modification time is used to "
                        "determine whether a job shall be executed, use the command "
                        "line option '--rerun-triggers mtime' (also see --help).\n"
                        "If you are sure that a change for a certain output file (say, <outfile>) won't "
                        "change the result (e.g. because you just changed the formatting of a script "
                        "or environment definition), you can also wipe its metadata to skip such a trigger via "
                        "'snakemake --cleanup-metadata <outfile>'. "
                    )
                    logger.info(
                        "Rules with provenance triggered jobs: "
                        + " ".join(jobs_to_rulenames(provenance_triggered_jobs))
                    )
                    logger.info("")

            has_checkpoint_jobs = any(self.dag.checkpoint_jobs)

            try:
                success = self.scheduler.schedule()
            except Exception as e:
                if self.dryrun:
                    log_provenance_info()
                raise e

            if (
                not self.remote_execution_settings.immediate_submit
                and not self.dryrun
                and self.exec_mode == ExecMode.DEFAULT
            ):
                self.dag.cleanup_workdir()

            if not dryrun_or_touch:
                async_run(self.dag.store_storage_outputs())
                async_run(self.dag.cleanup_storage_objects())

            if success:
                if self.dryrun:
                    if len(self.dag):
                        logger.run_info("\n".join(self.dag.stats()))
                        self.dag.print_reasons()
                        log_provenance_info()
                    logger.info("")
                    logger.info(
                        "This was a dry-run (flag -n). The order of jobs "
                        "does not reflect the order of execution."
                    )
                    if has_checkpoint_jobs:
                        logger.info(
                            "The run involves checkpoint jobs, "
                            "which will result in alteration of the DAG of "
                            "jobs (e.g. adding more jobs) after their completion."
                        )
                else:
                    logger.logfile_hint()
                if not self.dryrun and not self.execution_settings.no_hooks:
                    self._onsuccess(logger.get_logfile())
            else:
                if not self.dryrun and not self.execution_settings.no_hooks:
                    self._onerror(logger.get_logfile())
                logger.logfile_hint()
                raise WorkflowError("At least one job did not complete successfully.")

    @property
    def current_basedir(self):
        """Basedir of currently parsed Snakefile."""
        assert self.included_stack
        snakefile = self.included_stack[-1]
        basedir = snakefile.get_basedir()
        if isinstance(basedir, LocalSourceFile):
            return basedir.abspath()
        else:
            return basedir

    def source_path(self, rel_path):
        """Return path to source file from work dir derived from given path relative to snakefile"""
        # TODO download to disk (use source cache) in case of remote file
        import inspect

        frame = inspect.currentframe()
        assert frame is not None and frame.f_back is not None
        calling_file = frame.f_back.f_code.co_filename

        if (
            self.included_stack
            and calling_file == self.included_stack[-1].get_path_or_uri()
        ):
            # called from current snakefile, we can try to keep the original source
            # file annotation
            # This will only work if the method is evaluated during parsing mode.
            # Otherwise, the stack can be empty already.
            path = self.current_basedir.join(rel_path)
            orig_path = path.get_path_or_uri()
        else:
            # heuristically determine path
            calling_dir = os.path.dirname(calling_file)
            path = smart_join(calling_dir, rel_path)
            orig_path = path

        return sourcecache_entry(
            self.sourcecache.get_path(infer_source_file(path)), orig_path
        )

    @property
    def snakefile(self):
        import inspect

        frame = inspect.currentframe()
        assert not (frame is None or frame.f_back is None)
        return frame.f_back.f_code.co_filename

    def register_envvars(self, *envvars):
        """
        Register environment variables that shall be passed to jobs.
        If used multiple times, union is taken.
        """
        invalid_envvars = [
            envvar
            for envvar in envvars
            if re.match(r"^\w+$", envvar, flags=re.ASCII) is None
        ]
        if invalid_envvars:
            raise WorkflowError(
                f"Invalid environment variables requested: {', '.join(map(repr, invalid_envvars))}. "
                "Environment variable names may only contain alphanumeric characters and the underscore. "
            )
        undefined = set(var for var in envvars if var not in os.environ)
        if self.check_envvars and undefined:
            raise WorkflowError(
                "The following environment variables are requested by the workflow but undefined. "
                "Please make sure that they are correctly defined before running Snakemake:\n"
                "{}".format("\n".join(undefined))
            )
        self.envvars.update(envvars)

    def include(
        self,
        snakefile,
        overwrite_default_target=False,
        print_compilation=False,
    ):
        """
        Include a snakefile.
        """
        basedir = self.current_basedir if self.included_stack else None
        snakefile = infer_source_file(snakefile, basedir)

        if not self.modifier.allow_rule_overwrite and snakefile in self.included:
            logger.info(f"Multiple includes of {snakefile} ignored")
            return
        self.included.append(snakefile)
        self.included_stack.append(snakefile)

        default_target = self.default_target
        linemap: Dict[int, int] = dict()
        self.linemaps[snakefile.get_path_or_uri()] = linemap
        code, rulecount = parse(
            snakefile,
            self,
            rulecount=self._rulecount,
            linemap=linemap,
        )
        self._rulecount = rulecount

        if print_compilation:
            print(code)
            return

        snakefile_path_or_uri = snakefile.get_basedir().get_path_or_uri()
        if (
            isinstance(snakefile, LocalSourceFile)
            and snakefile_path_or_uri not in sys.path
        ):
            # insert the current directory into sys.path
            # this allows to import modules from the workflow directory
            sys.path.insert(0, snakefile.get_basedir().get_path_or_uri())

        exec(compile(code, snakefile.get_path_or_uri(), "exec"), self.globals)

        if not overwrite_default_target:
            self.default_target = default_target
        self.included_stack.pop()

    def onstart(self, func):
        """Register onstart function."""
        self._onstart = func

    def onsuccess(self, func):
        """Register onsuccess function."""
        self._onsuccess = func

    def onerror(self, func):
        """Register onerror function."""
        self._onerror = func

    def global_wildcard_constraints(self, **content):
        """Register global wildcard constraints."""
        self.modifier.wildcard_constraints.update(content)
        # update all rules so far
        for rule in self.modifier.rules:
            rule.update_wildcard_constraints()

    def scattergather(self, **content):
        """Register scattergather defaults."""
        self._scatter.update(content)
        self._scatter.update(self.resource_settings.overwrite_scatter)

        # add corresponding wildcard constraint
        self.global_wildcard_constraints(scatteritem=r"\d+-of-\d+")

        def func(key, *args, **wildcards):
            n = self._scatter[key]
            return expand(
                *args,
                scatteritem=map(f"{{}}-of-{n}".format, range(1, n + 1)),
                **wildcards,
            )

        for key in content:
            setattr(self.globals["scatter"], key, partial(func, key))
            setattr(self.globals["gather"], key, partial(func, key))

    def resourcescope(self, **content):
        """Register resource scope defaults"""
        self.resource_scopes.update(content)
        self.resource_scopes.update(self.resource_settings.overwrite_resource_scopes)

    def workdir(self, workdir):
        """Register workdir."""
        if self.overwrite_workdir is None:
            self._workdir_handler = WorkdirHandler(Path(workdir))
            self._workdir_handler.change_to()

    def configfile(self, fp):
        """Update the global config with data from the given file."""
        from snakemake.common.configfile import load_configfile

        if not self.modifier.skip_configfile:
            if os.path.exists(fp):
                self.configfiles.append(fp)
                c = load_configfile(fp)
                update_config(self.config, c)
                if self.config_settings.overwrite_config:
                    logger.info(
                        "Config file {} is extended by additional config specified via the command line.".format(
                            fp
                        )
                    )
                    update_config(self.config, self.config_settings.overwrite_config)
            elif not self.overwrite_configfiles:
                fp_full = os.path.abspath(fp)
                raise WorkflowError(
                    f"Workflow defines configfile {fp} but it is not present or accessible (full checked path: {fp_full})."
                )
            else:
                # CLI configfiles have been specified, do not throw an error but update with their values
                update_config(self.config, self.config_settings.overwrite_config)

    def set_pepfile(self, path):
        try:
            import peppy
        except ImportError:
            raise WorkflowError("For PEP support, please install peppy.")

        self.pepfile = str(path)
        self.globals["pep"] = peppy.Project(self.pepfile)

    def pepschema(self, schema):
        try:
            import eido
        except ImportError:
            raise WorkflowError("For PEP schema support, please install eido.")

        schema = str(schema)

        if is_local_file(schema) and not os.path.isabs(schema):
            # schema is relative to current Snakefile
            schema = self.current_basedir.join(schema).get_path_or_uri()
        if self.pepfile is None:
            raise WorkflowError("Please specify a PEP with the pepfile directive.")
        eido.validate_project(project=self.globals["pep"], schema=schema)

    def report(self, path):
        """Define a global report description in .rst format."""
        if not self.modifier.skip_global_report_caption:
            self.report_text = self.current_basedir.join(path)

    @property
    def config(self):
        return self.globals["config"]

    def ruleorder(self, *rulenames):
        self._ruleorder.add(*map(self.modifier.modify_rulename, rulenames))

    def localrules(self, *rulenames):
        self._localrules.update(rulenames)

    def rule(
        self, name: str | None = None, lineno=None, snakefile=None, checkpoint=False
    ):
        # choose a name for an unnamed rule
        orig_name = name or str(len(self._rules) + 1)

        if self.modifier.skip_rule(orig_name):

            def decorate(ruleinfo: RuleInfo):
                self.modifier.ruleinfos[orig_name] = ruleinfo
                ruleinfo.apply_modifier(self.modifier)
                # do nothing, ignore rule
                return ruleinfo.func

            return decorate

        # Optionally let the modifier change the rulename.
        rule = Rule(orig_name, self.modifier, lineno=lineno, snakefile=snakefile)
        # Register rule under its original name.
        # Modules using this snakefile as a module, will register it additionally under their
        # requested name.
        self.modifier.rule_proxies._register_rule(orig_name, RuleProxy(rule))

        rule.is_checkpoint = checkpoint
        rule.basedir = self.current_basedir
        rule.module_globals = self.modifier.globals

        self.add_rule(
            rule, checkpoint, allow_overwrite=self.modifier.allow_rule_overwrite
        )
        # handle default resources
        if self.resource_settings.default_resources is not None:
            rule.resources = copy.deepcopy(
                self.resource_settings.default_resources.parsed
            )
        else:
            rule.resources = dict()
        # Always require one node
        rule.resources["_nodes"] = 1

        if checkpoint:
            self.globals["checkpoints"].register(rule, fallback_name=orig_name)

        def decorate(ruleinfo: RuleInfo):  # type: ignore[no-redef]
            # If requested, modify ruleinfo via the modifier.
            ruleinfo.apply_modifier(self.modifier)

            if ruleinfo.name:
                del self._rules[rule.name]
                rule.name = ruleinfo.name
                self._rules[ruleinfo.name] = rule
            ruleinfo.func.__name__ = f"__{rule.name}"
            ruleinfo.update_rule(rule)
            ruleinfo.update_rule_settings(
                rule,
                resource_settings=self.resource_settings,
                workflow_settings=self.workflow_settings,
                container_img=self.global_container_img,
                is_containerized=self.global_is_containerized,
            )

            if ruleinfo.handover:
                if not ruleinfo.resources:
                    # give all available resources to the rule
                    rule.resources.update(
                        {
                            name: val
                            for name, val in self.global_resources.items()
                            if val is not None
                        }
                    )
                # This becomes a local rule, which might spawn jobs to a cluster,
                # depending on its configuration (e.g. nextflow config).
                self._localrules.add(rule.name)
                rule.is_handover = True

            if ruleinfo.cache and not (
                ruleinfo.cache is True
                or ruleinfo.cache == "omit-software"
                or ruleinfo.cache == "all"
            ):
                raise WorkflowError(
                    "Invalid value for cache directive. Use 'all' or 'omit-software'.",
                    rule=rule,
                )
            self.cache_rules[rule.name] = (
                "all" if ruleinfo.cache is True else ruleinfo.cache
            )

            if ruleinfo.default_target is True:
                self.default_target = rule.name
            elif ruleinfo.default_target is not False:
                raise WorkflowError(
                    "Invalid argument for 'default_target:' directive. Only True allowed. "
                    "Do not use the directive for rules that shall not be the default target. ",
                    rule=rule,
                )

            if ruleinfo.localrule is True:
                self._localrules.add(rule.name)

            self.globals[ruleinfo.func.__name__] = ruleinfo.func

            return ruleinfo.func

        return decorate

    @rule_decorate
    def docstring(self, ruleinfo: RuleInfo, string: str):
        ruleinfo.docstring = string.strip()

    def input(self, *paths, **kwpaths):
        def decorate(ruleinfo: RuleInfo):
            ruleinfo.input = InOutput(paths, kwpaths, self.modifier.path_modifier)
            return ruleinfo

        return decorate

    def output(self, *paths, **kwpaths):
        def decorate(ruleinfo: RuleInfo):
            ruleinfo.output = InOutput(paths, kwpaths, self.modifier.path_modifier)
            return ruleinfo

        return decorate

    def params(self, *params, **kwparams):

        def decorate(ruleinfo: RuleInfo):
            ruleinfo.params = (params, kwparams)
            return ruleinfo

        return decorate

    def register_wildcard_constraints(
        self, *wildcard_constraints, **kwwildcard_constraints
    ):
        def decorate(ruleinfo: RuleInfo):
            ruleinfo.wildcard_constraints = (
                wildcard_constraints,
                kwwildcard_constraints,
            )
            return ruleinfo

        return decorate

    @rule_decorate
    def cache_rule(self, ruleinfo: RuleInfo, cache: bool):
        ruleinfo.cache = cache

    @rule_decorate
    def default_target_rule(self, ruleinfo: RuleInfo, value):
        ruleinfo.default_target = value

    @rule_decorate
    def localrule(self, ruleinfo: RuleInfo, value):
        ruleinfo.localrule = value

    @rule_decorate
    def message(self, ruleinfo: RuleInfo, message):
        ruleinfo.message = message

    @rule_decorate
    def benchmark(self, ruleinfo: RuleInfo, benchmark):
        ruleinfo.benchmark = InOutput(benchmark, {}, self.modifier.path_modifier)

    @rule_decorate
    def conda(self, ruleinfo: RuleInfo, conda_env):
        ruleinfo.conda_env = conda_env

    def global_conda(self, conda_env):
        assert self.deployment_settings is not None
        if DeploymentMethod.CONDA in self.deployment_settings.deployment_method:
            from conda_inject import PackageManager, inject_env_file

            try:
                package_manager = PackageManager[
                    self.deployment_settings.conda_frontend.upper()
                ]
            except KeyError as err:
                raise WorkflowError(
                    f"Chosen conda frontend {self.deployment_settings.conda_frontend} is not supported by conda-inject."
                ) from err

            # Handle relative path
            if not isinstance(conda_env, SourceFile):
                if is_local_file(conda_env) and not os.path.isabs(conda_env):
                    # Conda env file paths are considered to be relative to the directory of the Snakefile
                    # hence we adjust the path accordingly.
                    # This is not necessary in case of receiving a SourceFile.
                    conda_env = self.current_basedir.join(conda_env)
                else:
                    # infer source file from unmodified uri or path
                    conda_env = infer_source_file(conda_env)

            logger.info(f"Injecting conda environment {conda_env.get_path_or_uri()}.")
            try:
                env = inject_env_file(
                    conda_env.get_path_or_uri(), package_manager=package_manager
                )
            except subprocess.CalledProcessError as e:
                raise WorkflowError(
                    f"Failed to inject conda environment {conda_env}: {e.stdout.decode()}",
                    e,
                )
            self.injected_conda_envs.append(env)

    @rule_decorate
    def container(self, ruleinfo: RuleInfo, container_img):
        # Explicitly set container_img to False if None is passed, indicating that
        # no container image shall be used, also not a global one.
        ruleinfo.container_img = container_img if container_img is not None else False
        ruleinfo.is_containerized = False

    @rule_decorate
    def containerized(self, ruleinfo: RuleInfo, container_img):
        ruleinfo.container_img = container_img
        ruleinfo.is_containerized = True

    @rule_decorate
    def envmodules(self, ruleinfo: RuleInfo, *env_modules):
        ruleinfo.env_modules = env_modules

    def global_container(self, container_img):
        self.global_container_img = container_img
        self.global_is_containerized = False

    def global_containerized(self, container_img):
        self.global_container_img = container_img
        self.global_is_containerized = True

    @rule_decorate
    def threads(self, ruleinfo: RuleInfo, threads):
        ruleinfo.threads = threads

    @rule_decorate
    def retries(self, ruleinfo: RuleInfo, retries):
        ruleinfo.retries = retries

    @rule_decorate
    def shadow(self, ruleinfo: RuleInfo, shadow_depth):
        ruleinfo.shadow_depth = shadow_depth

    def resources(self, *args, **resources):
        def decorate(ruleinfo: RuleInfo):
            ruleinfo.resources = (args, resources)
            return ruleinfo

        return decorate

    @rule_decorate
    def priority(self, ruleinfo: RuleInfo, priority):
        ruleinfo.priority = priority

    @rule_decorate
    def group(self, ruleinfo: RuleInfo, group):
        ruleinfo.group = group

    def log(self, *logs, **kwlogs):
        def decorate(ruleinfo):
            ruleinfo.log = InOutput(logs, kwlogs, self.modifier.path_modifier)
            return ruleinfo

        return decorate

    @rule_decorate
    def handover(self, ruleinfo: RuleInfo, value):
        ruleinfo.handover = value

    @rule_decorate
    def shellcmd(self, ruleinfo: RuleInfo, cmd):
        ruleinfo.shellcmd = cmd

    @rule_decorate
    def script(self, ruleinfo: RuleInfo, script):
        ruleinfo.script = script

    @rule_decorate
    def notebook(self, ruleinfo: RuleInfo, notebook):
        ruleinfo.notebook = notebook

    @rule_decorate
    def wrapper(self, ruleinfo: RuleInfo, wrapper):
        ruleinfo.wrapper = wrapper

    @rule_decorate
    def template_engine(self, ruleinfo: RuleInfo, template_engine):
        ruleinfo.template_engine = template_engine

    @rule_decorate
    def cwl(self, ruleinfo: RuleInfo, cwl):
        ruleinfo.cwl = cwl

    @rule_decorate
    def norun(self, ruleinfo: RuleInfo):
        ruleinfo.norun = True

    @rule_decorate
    def name(self, ruleinfo: RuleInfo, name):
        ruleinfo.name = name

    def run(self, func):
        return RuleInfo(func)

    def module(
        self,
        name,
        snakefile: str | None = None,
        meta_wrapper=None,
        config=None,
        skip_validation=False,
        replace_prefix=None,
        prefix=None,
    ):
        self.modules[name] = module = ModuleInfo(
            self,
            name,
            snakefile=snakefile,
            meta_wrapper=meta_wrapper,
            config=config,
            skip_validation=skip_validation,
            replace_prefix=replace_prefix,
            prefix=prefix,
        )
        module.load_ruleinfos(
            # do not overwrite existing report text via module
            skip_global_report_caption=self.report_text
            is not None,
        )

    def userule(
        self,
        rules: list[str],
        from_module: str | None = None,
        exclude_rules: list[str] | None = None,
        name_modifier: str | None = None,
        lineno: int | None = None,
    ):
        """
        Decorate a function or Ruleinfo object

        ```snakemake
        use rule xxx as yyy
        ```
        will be compiled as:
        ```python
        @workflow.userule(rules=["xxx"], name_modifier="yyy")
        def __userule_xxx___():
            pass
        ```

        comparatively
        ```snakemake
        use rule xxx as yyy with:
            params:
                a=1,
        ```
        will be compiled as:
        ```python
        @workflow.userule(rules=["xxx"], name_modifier="yyy")
        @workflow.params(a=1)
        @workflow.run
        def __userule_xxx___():
            pass
        ```
        where `workflow.run(lambda: None)` return a `RuleInfo` and not a callable
        """
        if from_module is not None:
            try:
                module = self.modules[from_module]
            except KeyError as err:
                raise WorkflowError(
                    f"Module {from_module} has not been registered with 'module' statement"
                    " before using it in 'use rule' statement."
                ) from err

            print(module.namespace._rules)

            def decorate(maybe_ruleinfo):
                ruleinfo = maybe_ruleinfo if not callable(maybe_ruleinfo) else None
                module.use_rules(
                    rules,
                    name_modifier,
                    exclude_rules=exclude_rules,
                    ruleinfo=ruleinfo,
                    # do not overwrite existing report text via module
                    skip_global_report_caption=self.report_text is not None,
                )

        else:
            if len(rules) > 1:
                raise WorkflowError(
                    "'use rule' statement from rule in the same module "
                    "must declare a single rule but multiple rules are declared."
                )
            # The parent use rule statement is specific for a different particular rule
            # else, we just add the rule into ruleinfo and skip it.
            if self.modifier.skip_rule(name_modifier):
                rule_whitelist: list | None = []
            else:
                rule_whitelist = None
            import copy

            orig_ruleinfo = copy.copy(self.modifier.get_ruleinfo(rules[0]))

            def decorate(maybe_ruleinfo):
                # local inheritance
                if rule_whitelist is not None and name_modifier is None:
                    # assert rule_whitelist == [], "means this rule is skipped"
                    # name_modifier is None means the rule is unnamed and will never be
                    # refered as rules.xxx anymore.
                    # Hense it can be ignored safely.
                    return
                ruleinfo = maybe_ruleinfo if not callable(maybe_ruleinfo) else None
                ruleinfos = self.modifier.ruleinfos
                with WorkflowModifier(
                    self,
                    parent_modifier=self.modifier,
                    resolved_rulename_modifier=get_name_modifier_func(
                        rules, name_modifier, parent_modifier=self.modifier
                    ),
                    rule_whitelist=rule_whitelist,
                    ruleinfo_overwrite=ruleinfo,
                ):
                    # A copy is necessary to avoid leaking modifications in case of multiple inheritance statements.
                    self.rule(
                        name=name_modifier,
                        lineno=lineno,
                        snakefile=self.included_stack[-1],
                    )(orig_ruleinfo)
                    ruleinfos |= self.modifier.ruleinfos

        return decorate

    @staticmethod
    def _empty_decorator(f):
        return f
