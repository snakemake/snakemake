from abc import ABC
from dataclasses import dataclass
from enum import Enum
import importlib
from pathlib import Path
from typing import Dict, List, Optional, Self, Set, Type
from snakemake import notebook
from snakemake.rules import Rule

from snakemake_interface_executor_plugins.utils import ExecMode
from snakemake_interface_executor_plugins.settings import (
    RemoteExecutionSettingsExecutorInterface,
    ConfigSettingsExecutorInterface,
    DeploymentSettingsExecutorInterface,
    ExecutionSettingsExecutorInterface,
    OutputSettingsExecutorInterface,
    ResourceSettingsExecutorInterface,
    StorageSettingsExecutorInterface,
    SettingsEnumBase,
    DeploymentMethod,
)

from snakemake.common import RERUN_TRIGGERS, RerunTrigger, dict_to_key_value_args
from snakemake.dag import Batch
from snakemake.exceptions import ApiError
from snakemake.io import load_configfile
from snakemake.resources import DefaultResources
from snakemake.utils import update_config
from snakemake.exceptions import WorkflowError
    

class RerunTrigger(SettingsEnumBase):
    MTIME = 0
    PARAMS = 1
    INPUT = 2
    SOFTWARE_ENV = 3
    CODE = 4


class SettingsBase(ABC):
    def __post_init__(self):
        self._check()

    def _check(self):
        pass

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
    workdir: Optional[Path] = None
    cache: Optional[List[str]] = None
    keep_going: bool = False
    debug: bool = False
    standalone: bool = False
    ignore_ambiguity: bool = False
    lock: bool = True
    ignore_incomplete: bool = False
    latency_wait: 3
    wait_for_files: Optional[List[str]] = None
    notemp: bool = False
    all_temp: bool = False
    keep_target_files: bool = False
    no_hooks: bool = False
    overwrite_shellcmd: Optional[str] = None
    retries: int = 0
    attempt: int = 1
    use_threads: bool = False
    shadow_prefix: Optional[Path] = None
    mode: ExecMode = ExecMode.default
    wrapper_prefix: Optional[str] = None
    keep_incomplete: bool = False
    keep_metadata: bool = True
    max_inventory_wait_time: int = 20
    edit_notebook: Optional[notebook.EditMode] = None
    cleanup_scripts: bool = True
    cleanup_metadata: List[Path] = []

@dataclass
class DAGSettings(SettingsBase):
    targets: Optional[List[str]] = None
    target_jobs: Set(str) = set()
    batch: Batch = None
    forcetargets: bool = False
    forceall: bool = False
    forcerun: List[str] = []
    until: List[str] = []
    omit_from: List[str] = []
    force_incomplete: bool = False
    allowed_rules: Set[str] = {}
    rerun_triggers: Set[RerunTrigger]=RerunTrigger.all()

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
                raise ApiError(f"Unknown default remote provider {self.default_remote_provider}.")
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
    deployment_method: Set[DeploymentMethod] = set()
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
    prioritytargets: Optional[List[str]] = None
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
            raise ApiError(
                "greediness must be >0 and <=1"
            )

@dataclass
class ResourceSettings(SettingsBase, ResourceSettingsExecutorInterface):
    cores: Optional[int] = None
    nodes: Optional[int] = None
    local_cores: Optional[int] = None
    max_threads: Optional[int] = None
    resources: Dict[str, int] = dict()
    overwrite_threads: Dict[str, int] = dict()
    overwrite_scatter: Dict[str, int] = dict()
    overwrite_resource_scopes: Dict[str, str] = dict()
    overwrite_resources: Dict[str, Dict[str, int]] = dict()
    default_resources: DefaultResources = DefaultResources(mode="bare")

    def _check(self):
        if self.cores is None and self.nodes is None:
            raise ApiError("You need to specify either --cores or --jobs.")

@dataclass
class ConfigSettings(SettingsBase, ConfigSettingsExecutorInterface):
    config: Dict[str, str] = dict()
    configfiles: Optional[List[Path]] = []
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
class OutputSettings(SettingsBase, OutputSettingsExecutorInterface):
    printshellcmds: bool = False
    nocolor: bool = False
    quiet: Optional[Quietness] = None
    debug_dag: bool = False
    verbose: bool = False
    show_failed_logs: bool = False
    log_handlers: List[object] = None
    keep_logger: bool = False


@dataclass
class PreemptibleRules:
    rules: Optional[Set[str]] = None
    all: bool = False

    def __post_init__(self):
        assert self.rules or self.all, "bug: either rules or all have to be set"
    
    def is_preemptible(self, rulename: str):
        return self.all or rulename in self.rules


@dataclass
class RemoteExecutionSettings(SettingsBase, RemoteExecutionSettingsExecutorInterface):
    jobname: str = "snakejob.{rulename}.{jobid}.sh"
    jobscript: Optional[Path] = None
    max_status_checks_per_second: int = 100
    container_image: Optional[str] = None
    preemptible_retries: Optional[int] = None
    preemptible_rules: Optional[PreemptibleRules] = None
    envvars: Optional[List[str]] = None
    immediate_submit: bool = False

@dataclass
class GroupSettings(SettingsBase):
    overwrite_groups: Dict[str, str] = dict()
    group_components: Dict[str, int] = dict()
    local_groupid: str = "local"