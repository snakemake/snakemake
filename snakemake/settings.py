from abc import ABC
from dataclasses import dataclass
from enum import Enum
import importlib
from pathlib import Path
from typing import Dict, List, Optional, Set

from snakemake_interface_executor_plugins.utils import ExecMode
from snakemake_interface_executor_plugins.settings import (
    RemoteExecutionSettingsExecutorInterface,
    ConfigSettingsExecutorInterface,
    DeploymentSettingsExecutorInterface,
    ExecutionSettingsExecutorInterface,
    OutputSettingsExecutorInterface,
    ResourceSettingsExecutorInterface,
    StorageSettingsExecutorInterface,
)

from snakemake.common import RERUN_TRIGGERS, dict_to_key_value_args
from snakemake.dag import Batch
from snakemake.exceptions import ApiError
from snakemake.io import load_configfile
from snakemake.resources import DefaultResources
from snakemake.utils import update_config


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
    keep_remote_local: bool = False
    keep_target_files: bool = False
    no_hooks: bool = False
    overwrite_shellcmd: Optional[str] = None
    updated_files: List[str] = []
    restart_times: int = 0
    attempt: int = 1
    use_threads: bool = False
    shadow_prefix: Optional[Path] = None
    mode: ExecMode = ExecMode.default
    wrapper_prefix: Optional[str] = None
    keep_incomplete: bool = False
    keep_metadata: bool = True
    max_inventory_wait_time: int = 20
    edit_notebook: Optional[Path] = None
    cleanup_scripts: bool = True
    cleanup_metadata: List[Path] = []

    def _check(self):
        # TODO move into API as immediate_submit has been moved to remote_execution_settings
        assert not self.immediate_submit or (
            self.immediate_submit and self.notemp
        ), "immediate_submit has to be combined with notemp (it does not support temp file handling)"

@dataclass
class DAGSettings(SettingsBase):
    targets: Optional[List[str]] = None
    target_jobs: Set(str) = set()
    batch: Batch = None
    forcetargets: bool = False
    forceall: bool = False
    forcerun: Optional[List[str]] = None
    until: Optional[List[str]] = None
    omit_from: Optional[List[str]] = None
    force_incomplete: bool = False
    allowed_rules: Optional[Set[str]] = None
    rerun_triggers: List[str]=RERUN_TRIGGERS

@dataclass
class StorageSettings(SettingsBase, StorageSettingsExecutorInterface):
    default_remote_provider: Optional[str] = None
    default_remote_prefix: Optional[str] = None
    assume_shared_fs: bool = True
    keep_remote_local: bool = False

    def __post_init__(self):
        self.default_remote_provider = self._get_default_remote_provider()

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


class CondaCleanupPkgs(Enum):
    TARBALLS = "tarballs"
    CACHE = "cache"

@dataclass
class DeploymentSettings(SettingsBase, DeploymentSettingsExecutorInterface):
    """
    Parameters
    ----------

    use_conda:
        use conda environments for each job (defined with conda directive of rules)
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
    use_conda: bool = False
    conda_prefix: Optional[Path] = None
    conda_cleanup_pkgs: Optional[CondaCleanupPkgs] = None
    conda_base_path: Optional[Path] = None
    conda_frontend: str = "mamba"
    conda_not_block_search_path_envvars: bool = False
    use_singularity: bool = False
    singularity_args: str = ""
    singularity_prefix: Optional[Path] = None
    use_env_modules: bool = False



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
    cores: Optional[int]
    nodes: Optional[int] = None
    local_cores: Optional[int] = None
    max_threads: Optional[int] = None
    resources: Dict[str, int] = dict()
    overwrite_threads: Dict[str, int] = dict()
    overwrite_scatter: Dict[str, int] = dict()
    overwrite_resource_scopes: Dict[str, str] = dict()
    overwrite_resources: Dict[str, Dict[str, int]] = dict()
    default_resources: DefaultResources = DefaultResources(mode="bare")

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

@dataclass
class OutputSettings(SettingsBase, OutputSettingsExecutorInterface):
    printshellcmds: bool = False
    nocolor: bool = False
    quiet: bool = False
    debug_dag: bool = False
    verbose: bool = False
    show_failed_logs: bool = False
    log_handler: List[object] = None
    keep_logger: bool = False


@dataclass
class RemoteExecutionSettings(SettingsBase, RemoteExecutionSettingsExecutorInterface):
    jobname: str = "snakejob.{rulename}.{jobid}.sh"
    jobscript: Optional[Path] = None
    max_status_checks_per_second: int = 100
    container_image: Optional[str] = None
    preemption_default: Optional[int] = None
    preemptible_rules: Optional[List[str]] = None
    envvars: Optional[List[str]] = None
    immediate_submit: bool = False

@dataclass
class GroupSettings(SettingsBase):
    overwrite_groups: Dict[str, str] = dict()
    group_components: Dict[str, int] = dict()
    local_groupid: str = "local"