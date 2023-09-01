from abc import ABC
from dataclasses import dataclass, field
from enum import Enum
import importlib
from pathlib import Path
from types import MappingProxyType
from typing import Optional
from collections.abc import Mapping, Sequence, Set

from snakemake_interface_common.exceptions import ApiError
from snakemake_interface_executor_plugins.settings import (
    RemoteExecutionSettingsExecutorInterface,
    DeploymentSettingsExecutorInterface,
    ExecutionSettingsExecutorInterface,
    StorageSettingsExecutorInterface,
    DeploymentMethod,
    ExecMode,
)
from snakemake_interface_common.settings import SettingsEnumBase

from snakemake.common import dict_to_key_value_args
from snakemake.common.configfile import load_configfile
from snakemake.resources import DefaultResources
from snakemake.utils import update_config
from snakemake.exceptions import WorkflowError


frozendict = lambda: MappingProxyType(dict())


class frozendict(dict):
    def __init__(self):
        super().__init__(dict())

    def __hash__(self):
        return hash(tuple(sorted(self.items())))


class RerunTrigger(SettingsEnumBase):
    MTIME = 0
    PARAMS = 1
    INPUT = 2
    SOFTWARE_ENV = 3
    CODE = 4


class ChangeType(SettingsEnumBase):
    CODE = 0
    INPUT = 1
    PARAMS = 2


class SettingsBase(ABC):
    def __post_init__(self):
        self._check()

    def _check(self):
        pass


class NotebookEditMode:
    def __init__(self, server_addr: Optional[str] = None, draft_only: bool = False):
        if server_addr is not None:
            self.ip, self.port = server_addr.split(":")
        self.draft_only = draft_only


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
    cache: Optional[Sequence[str]] = None
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
    mode: ExecMode = ExecMode.DEFAULT
    wrapper_prefix: Optional[str] = None
    keep_incomplete: bool = False
    keep_metadata: bool = True
    edit_notebook: Optional[NotebookEditMode] = None
    cleanup_scripts: bool = True
    cleanup_metadata: Sequence[Path] = tuple()


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


@dataclass
class DAGSettings(SettingsBase):
    targets: Set[str] = frozenset()
    target_jobs: Set[str] = frozenset()
    target_files_omit_workdir_adjustment: bool = False
    batch: Batch = None
    forcetargets: bool = False
    forceall: bool = False
    forcerun: Set[str] = frozenset()
    until: Set[str] = frozenset()
    omit_from: Set[str] = frozenset()
    force_incomplete: bool = False
    allowed_rules: Set[str] = frozenset()
    rerun_triggers: Set[RerunTrigger] = RerunTrigger.all()
    max_inventory_wait_time: int = 20

    def _check(self):
        if self.batch is not None and self.forceall:
            raise WorkflowError(
                "--batch may not be combined with --forceall, because recomputed upstream "
                "jobs in subsequent batches may render already obtained results outdated."
            )


@dataclass
class StorageSettings(SettingsBase, StorageSettingsExecutorInterface):
    default_remote_provider: Optional[str] = None
    default_remote_prefix: Optional[str] = None
    assume_shared_fs: bool = True
    keep_remote_local: bool = False
    notemp: bool = False
    all_temp: bool = False

    def __post_init__(self):
        self.default_remote_provider = self._get_default_remote_provider()
        super().__post_init__()

    def _get_default_remote_provider(self):
        if self.default_remote_provider is not None:
            try:
                rmt = importlib.import_module(
                    "snakemake.remote." + self.default_remote_provider
                )
            except ImportError as e:
                raise ApiError(
                    f"Unknown default remote provider {self.default_remote_provider}."
                )
            if rmt.RemoteProvider.supports_default:
                return rmt.RemoteProvider(
                    keep_local=self.keep_remote_local, is_default=True
                )
            else:
                raise ApiError(
                    "Remote provider {} does not (yet) support to "
                    "be used as default provider."
                )


class CondaCleanupPkgs(SettingsEnumBase):
    TARBALLS = 0
    CACHE = 1


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

    deployment_method: Set[DeploymentMethod] = frozenset()
    conda_prefix: Optional[Path] = None
    conda_cleanup_pkgs: Optional[CondaCleanupPkgs] = None
    conda_base_path: Optional[Path] = None
    conda_frontend: str = "mamba"
    conda_not_block_search_path_envvars: bool = False
    apptainer_args: str = ""
    apptainer_prefix: Optional[Path] = None


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
        set the greediness of scheduling. This value between 0 and 1 determines how careful jobs are selected for execution. The default value (0.5 if prioritytargets are used, 1.0 else) provides the best speed and still acceptable scheduling quality.
    """

    prioritytargets: Set[str] = frozenset()
    scheduler: str = "ilp"
    ilp_solver: Optional[str] = None
    solver_path: Optional[Path] = None
    greediness: Optional[float] = None
    max_jobs_per_second: int = 10

    def __post_init__(self):
        self.greediness = self._get_greediness()

    def _get_greediness(self):
        if self.greediness is None:
            return 0.5 if self.prioritytargets else 1.0
        else:
            return self.greediness

    def _check(self):
        if not (0 < self.greedyness <= 1.0):
            raise ApiError("greediness must be >0 and <=1")


@dataclass
class ResourceSettings(SettingsBase):
    cores: Optional[int] = None
    nodes: Optional[int] = None
    local_cores: Optional[int] = None
    max_threads: Optional[int] = None
    resources: Mapping[str, int] = frozendict()
    overwrite_threads: Mapping[str, int] = frozendict()
    overwrite_scatter: Mapping[str, int] = frozendict()
    overwrite_resource_scopes: Mapping[str, str] = frozendict()
    overwrite_resources: Mapping[str, Mapping[str, int]] = frozendict()
    default_resources: Optional[DefaultResources] = None


@dataclass
class ConfigSettings(SettingsBase):
    config: Mapping[str, str] = frozendict()
    configfiles: Sequence[Path] = tuple()
    config_args: Optional[str] = None

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

    def _get_configfiles(self):
        return list(map(Path.absolute, self.configfiles))

    def _get_config_args(self):
        if self.config_args is None:
            return dict_to_key_value_args(self.config)
        else:
            return self.config_args


class Quietness(SettingsEnumBase):
    RULES = 0
    PROGRESS = 1
    ALL = 2


@dataclass
class OutputSettings(SettingsBase):
    printshellcmds: bool = False
    nocolor: bool = False
    quiet: Optional[Set[Quietness]] = None
    debug_dag: bool = False
    verbose: bool = False
    show_failed_logs: bool = False
    log_handlers: Sequence[object] = tuple()
    keep_logger: bool = False


@dataclass
class PreemptibleRules:
    rules: Set[str] = frozenset()
    all: bool = False

    def is_preemptible(self, rulename: str):
        return self.all or rulename in self.rules


@dataclass
class RemoteExecutionSettings(SettingsBase, RemoteExecutionSettingsExecutorInterface):
    jobname: str = "snakejob.{rulename}.{jobid}.sh"
    jobscript: Optional[Path] = None
    max_status_checks_per_second: int = 100
    container_image: Optional[str] = None
    preemptible_retries: Optional[int] = None
    preemptible_rules: PreemptibleRules = field(default_factory=PreemptibleRules)
    envvars: Sequence[str] = tuple()
    immediate_submit: bool = False


@dataclass
class GroupSettings(SettingsBase):
    overwrite_groups: Mapping[str, str] = frozendict()
    group_components: Mapping[str, int] = frozendict()
    local_groupid: str = "local"
