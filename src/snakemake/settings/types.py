from abc import ABC
from dataclasses import dataclass, field
import os
from pathlib import Path
import re
from typing import Any, Optional, Union
from typing import Mapping, Sequence, Set

import immutables

from snakemake.common.typing import AnySet
from snakemake_interface_common.exceptions import ApiError
from snakemake_interface_executor_plugins.settings import (
    RemoteExecutionSettingsExecutorInterface,
    DeploymentSettingsExecutorInterface,
    ExecutionSettingsExecutorInterface,
    StorageSettingsExecutorInterface,
    DeploymentMethod,
    ExecMode,
    SharedFSUsage,
)
from snakemake_interface_logger_plugins.settings import (
    LogHandlerSettingsBase,
    OutputSettingsLoggerInterface,
)

from snakemake.common import (
    dict_to_key_value_args,
    expand_vars_and_user,
    get_container_image,
)
from snakemake.common.configfile import load_configfile
from snakemake.resources import DefaultResources
from snakemake.utils import update_config
from snakemake.exceptions import WorkflowError
from snakemake.settings.enums import (
    RerunTrigger,
    ChangeType,
    CondaCleanupPkgs,
    Quietness,
    StrictDagEvaluation,
    PrintDag,
)


class SettingsBase(ABC):
    def __post_init__(self):
        self._check()

    def _check(self):
        # by default, nothing to check
        # override this method in subclasses if needed
        pass


class NotebookEditMode:
    def __init__(self, server_addr: Optional[str] = None, draft_only: bool = False):
        if server_addr is not None:
            self.ip, self.port = server_addr.split(":")
        self.draft_only = draft_only


class MaxJobsPerTimespan:
    arg_re = re.compile(r"(?P<count>\d+)/(?P<timespan>\d+(h|m|s|ms|w|d))")

    def __init__(self, max_jobs: int, timespan: Optional[str] = None):
        import humanfriendly

        self.max_jobs = max_jobs
        if timespan is not None:
            self.timespan = humanfriendly.parse_timespan(timespan)
        else:
            self.timespan = 1

    @classmethod
    def parse_choice(cls, arg: str):
        m = cls.arg_re.match(arg)
        if m is None:
            raise WorkflowError(
                "Invalid max jobs per timespan definition. "
                "Must be of the form <max_jobs>/<timespan> with <max_jobs> being an "
                "integer, and <timespan> being an integer with "
                f"unit h, m, s ms, w, d. Given instead: {arg}"
            )
        max_jobs, timespan = m.group("count"), m.group("timespan")
        max_jobs = int(max_jobs)
        return cls(max_jobs, timespan=timespan)


@dataclass
class ExecutionSettings(SettingsBase, ExecutionSettingsExecutorInterface):
    """
    Parameters
    ----------

    batch:
        whether to compute only a partial DAG, defined by the given Batch object
    cache:
        list of rules to cache
    cores:
        the number of provided cores (ignored when using cluster/cloud support)
    nodes:
        the number of provided cluster nodes (ignored without cluster/cloud support)
    local_cores:
        the number of provided local cores if in cluster mode (ignored without cluster/cloud support)
    """

    latency_wait: int = 3
    keep_going: bool = False
    debug: bool = False
    standalone: bool = False
    ignore_ambiguity: bool = False
    lock: bool = True
    ignore_incomplete: bool = False
    wait_for_files: Sequence[str] = tuple()
    no_hooks: bool = False
    retries: int = 0
    attempt: int = 1
    use_threads: bool = False
    shadow_prefix: Optional[Path] = None
    keep_incomplete: bool = False
    keep_metadata: bool = True
    edit_notebook: Optional[NotebookEditMode] = None
    cleanup_scripts: bool = True
    queue_input_wait_time: int = 10


@dataclass
class WorkflowSettings(SettingsBase):
    wrapper_prefix: Optional[str] = None
    exec_mode: ExecMode = ExecMode.DEFAULT
    cache: Optional[Sequence[str]] = None
    consider_ancient: Mapping[str, AnySet[Union[str, int]]] = field(
        default_factory=dict
    )


class Batch:
    """Definition of a batch for calculating only a partial DAG."""

    def __init__(self, rulename: str, idx: int, batches: int):
        assert idx <= batches
        assert idx > 0
        self.rulename = rulename
        self.idx = idx
        self.batches = batches

    def get_batch(self, items: list):
        """Return the defined batch of the given items.
        Items are usually input files."""
        # make sure that we always consider items in the same order
        if len(items) < self.batches:
            raise WorkflowError(
                "Batching rule {} has less input files than batches. "
                "Please choose a smaller number of batches.".format(self.rulename)
            )
        items = sorted(items)

        # we can equally split items using divmod:
        # len(items) = (self.batches * quotient) + remainder
        # Because remainder always < divisor (self.batches),
        # each batch will be equal to quotient + (1 or 0 item)
        # from the remainder
        k, m = divmod(len(items), self.batches)

        # self.batch is one-based, hence we have to subtract 1
        idx = self.idx - 1

        # First n batches will have k (quotient) items +
        # one item from the remainder (m). Once we consume all items
        # from the remainder, last batches only contain k items.
        i = idx * k + min(idx, m)
        batch_len = (idx + 1) * k + min(idx + 1, m)

        if self.is_final:
            # extend the last batch to cover rest of list
            return items[i:]
        else:
            return items[i:batch_len]

    @property
    def is_final(self):
        return self.idx == self.batches

    def __str__(self):
        return f"{self.idx}/{self.batches} (rule {self.rulename})"

    def __eq__(self, other):
        return (
            self.rulename == other.rulename
            and self.idx == other.idx
            and self.batches == other.batches
        )


@dataclass
class DAGSettings(SettingsBase):
    targets: AnySet[str] = frozenset()
    target_jobs: AnySet[str] = frozenset()
    target_files_omit_workdir_adjustment: bool = False
    batch: Optional[Batch] = None
    forcetargets: bool = False
    forceall: bool = False
    forcerun: AnySet[str] = frozenset()
    until: AnySet[str] = frozenset()
    omit_from: AnySet[str] = frozenset()
    force_incomplete: bool = False
    allowed_rules: AnySet[str] = frozenset()
    rerun_triggers: AnySet[RerunTrigger] = RerunTrigger.all()
    max_inventory_wait_time: int = 20
    trust_io_cache: bool = False
    max_checksum_file_size: int = 1000000
    strict_evaluation: AnySet[StrictDagEvaluation] = frozenset()
    print_dag_as: PrintDag = PrintDag.DOT
    # strict_functions_evaluation: bool = False
    # strict_cycle_evaluation: bool = False
    # strict_wildcards_recursion_evaluation: bool = False

    def _check(self):
        if self.batch is not None and self.forceall:
            raise WorkflowError(
                "--batch may not be combined with --forceall, because recomputed upstream "
                "jobs in subsequent batches may render already obtained results outdated."
            )


@dataclass
class StorageSettings(SettingsBase, StorageSettingsExecutorInterface):
    default_storage_provider: Optional[str] = None
    default_storage_prefix: Optional[str] = None
    shared_fs_usage: AnySet[SharedFSUsage] = SharedFSUsage.all()
    keep_storage_local: bool = False
    retrieve_storage: bool = True
    local_storage_prefix: Path = Path(".snakemake/storage")
    remote_job_local_storage_prefix: Optional[Path] = None
    notemp: bool = False
    all_temp: bool = False
    unneeded_temp_files: AnySet[str] = frozenset()
    wait_for_free_local_storage: Optional[int] = None

    def __post_init__(self):
        if self.remote_job_local_storage_prefix is None:
            self.remote_job_local_storage_prefix = self.local_storage_prefix


@dataclass
class DeploymentSettings(SettingsBase, DeploymentSettingsExecutorInterface):
    """
    Parameters
    ----------

    deployment_method
        deployment method to use (CONDA, APPTAINER, ENV_MODULES)
    conda_prefix:
        the directory in which conda environments will be created (default None)
    conda_cleanup_pkgs:
        whether to clean up conda tarballs after env creation (default None), valid values: "tarballs", "cache"
    conda_create_envs_only:
        if specified, only builds the conda environments specified for each job, then exits.
    list_conda_envs:
        list conda environments and their location on disk.
    conda_base_path:
        Path to conda base environment (this can be used to overwrite the search path for conda, mamba, and activate).
    """

    deployment_method: AnySet[DeploymentMethod] = frozenset()
    conda_prefix: Optional[Path] = None
    conda_cleanup_pkgs: Optional[CondaCleanupPkgs] = None
    conda_base_path: Optional[Path] = None
    conda_frontend: str = "conda"
    conda_not_block_search_path_envvars: bool = False
    apptainer_args: str = ""
    apptainer_prefix: Optional[Path] = None

    def imply_deployment_method(self, method: DeploymentMethod):
        self.deployment_method = set(self.deployment_method)
        self.deployment_method.add(method)

    def __post_init__(self):
        from snakemake.logging import logger

        if self.apptainer_prefix is None:
            self.apptainer_prefix = os.environ.get("APPTAINER_CACHEDIR", None)
        self.apptainer_prefix = expand_vars_and_user(self.apptainer_prefix)
        self.conda_prefix = expand_vars_and_user(self.conda_prefix)
        if self.conda_frontend != "conda":
            logger.warning(
                "Support for alternative conda frontends has been deprecated in "
                "favor of simpler support and code base. "
                "This should not cause issues since current conda releases rely on "
                "fast solving via libmamba. "
                f"Ignoring the alternative conda frontend setting ({self.conda_frontend})."
            )
            self.conda_frontend = "conda"


@dataclass
class SchedulingSettings(SettingsBase):
    """
    Parameters
    ----------

    prioritytargets:
        list of targets that shall be run with maximum priority (default [])
    scheduler:
        Select scheduling algorithm (default ilp, allowed: ilp, greedy)
    ilp_solver:
        Set solver for ilp scheduler.
    greediness:
        Set the greediness of scheduling. This value, between 0 and 1, determines how careful jobs are selected for execution. The default value (0.5 if prioritytargets are used, 1.0 else) provides the best speed and still acceptable scheduling quality.
    subsample:
        Set the number of jobs to be considered for scheduling. If number of ready jobs is greater than this value, this number of jobs is randomly chosen for scheduling; if number of ready jobs is lower, this option has no effect. This can be useful on very large DAGs, where the scheduler can take some time selecting which jobs to run."
    """

    prioritytargets: AnySet[str] = frozenset()
    scheduler: str = "ilp"
    ilp_solver: Optional[str] = None
    solver_path: Optional[Path] = None
    greediness: Optional[float] = None
    subsample: Optional[int] = None
    max_jobs_per_second: Optional[int] = None
    max_jobs_per_timespan: Optional[MaxJobsPerTimespan] = None

    def __post_init__(self):
        self.greediness = self._get_greediness()
        if self.max_jobs_per_second is not None and self.max_jobs_per_timespan is None:
            self.max_jobs_per_timespan = MaxJobsPerTimespan(self.max_jobs_per_second)

    def _get_greediness(self):
        if self.greediness is None:
            return 0.5 if self.prioritytargets else 1.0
        else:
            return self.greediness

    def _check(self):
        if not (0 <= self.greediness <= 1.0):
            raise ApiError("greediness must be >=0 and <=1")
        if self.subsample:
            if not isinstance(self.subsample, int) or self.subsample < 1:
                raise ApiError("subsample must be a positive integer")


@dataclass
class ResourceSettings(SettingsBase):
    cores: Optional[int] = None
    nodes: Optional[int] = None
    local_cores: Optional[int] = None
    max_threads: Optional[int] = None
    resources: Mapping[str, int] = immutables.Map()
    overwrite_threads: Mapping[str, int] = immutables.Map()
    overwrite_scatter: Mapping[str, int] = immutables.Map()
    overwrite_resource_scopes: Mapping[str, str] = immutables.Map()
    overwrite_resources: Mapping[str, Mapping[str, Any]] = immutables.Map()
    default_resources: Optional[DefaultResources] = None

    def __post_init__(self):
        if self.default_resources is None:
            self.default_resources = DefaultResources(mode="bare")


@dataclass
class ConfigSettings(SettingsBase):
    config: Mapping[str, str] = immutables.Map()
    configfiles: Sequence[Path] = tuple()
    config_args: Optional[str] = None
    replace_workflow_config: bool = False

    def __post_init__(self):
        self.overwrite_config = self._get_overwrite_config()
        self.configfiles = self._get_configfiles()
        self.config_args = self._get_config_args()

    def _get_overwrite_config(self):
        overwrite_config = dict()
        if self.configfiles:
            for f in self.configfiles:
                update_config(overwrite_config, load_configfile(f))
        if self.config:
            update_config(overwrite_config, self.config)
        return overwrite_config

    def _get_configfiles(self):
        return list(map(Path.absolute, self.configfiles))

    def _get_config_args(self):
        if self.config_args is None:
            return dict_to_key_value_args(self.config, repr_obj=True)
        else:
            return self.config_args


@dataclass
class OutputSettings(SettingsBase, OutputSettingsLoggerInterface):
    dryrun: bool = False
    printshellcmds: bool = False
    nocolor: bool = False
    quiet: Optional[AnySet[Quietness]] = None
    debug_dag: bool = False
    verbose: bool = False
    show_failed_logs: bool = False
    log_handler_settings: Mapping[str, LogHandlerSettingsBase] = immutables.Map()
    keep_logger: bool = False
    stdout: bool = False
    benchmark_extended: bool = False


@dataclass
class PreemptibleRules:
    rules: AnySet[str] = frozenset()
    all: bool = False

    def is_preemptible(self, rulename: str):
        return self.all or rulename in self.rules


@dataclass
class RemoteExecutionSettings(SettingsBase, RemoteExecutionSettingsExecutorInterface):
    jobname: str = "snakejob.{rulename}.{jobid}.sh"
    jobscript: Optional[Path] = None
    max_status_checks_per_second: float = 100.0
    seconds_between_status_checks: int = 10
    container_image: str = get_container_image()
    preemptible_retries: Optional[int] = None
    preemptible_rules: PreemptibleRules = field(default_factory=PreemptibleRules)
    envvars: Sequence[str] = tuple()
    immediate_submit: bool = False
    precommand: Optional[str] = None
    job_deploy_sources: bool = True


@dataclass
class GroupSettings(SettingsBase):
    overwrite_groups: Mapping[str, str] = immutables.Map()
    group_components: Mapping[str, int] = immutables.Map()
    local_groupid: str = "local"
