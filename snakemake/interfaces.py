__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import sys

from abc import ABC, abstractmethod
from typing import Optional, List, Dict


class ExecutorJobInterface(ABC):
    @property
    @abstractmethod
    def name(self):
        ...

    @property
    @abstractmethod
    def jobid(self):
        ...

    @abstractmethod
    def logfile_suggestion(self, prefix: str) -> str:
        ...

    @abstractmethod
    def is_group(self):
        ...

    @abstractmethod
    def log_info(self, skip_dynamic=False):
        ...

    @abstractmethod
    def log_error(self, msg=None, **kwargs):
        ...

    @abstractmethod
    def remove_existing_output(self):
        ...

    @abstractmethod
    def download_remote_input(self):
        ...

    @abstractmethod
    def properties(self, omit_resources=["_cores", "_nodes"], **aux_properties):
        ...

    @property
    @abstractmethod
    def resources(self):
        ...

    @abstractmethod
    def check_protected_output(self):
        ...

    @property
    @abstractmethod
    def is_local(self):
        ...

    @property
    @abstractmethod
    def is_branched(self):
        ...

    @property
    @abstractmethod
    def is_updated(self):
        ...

    @property
    @abstractmethod
    def output(self):
        ...

    @abstractmethod
    def register(self):
        ...

    @abstractmethod
    def postprocess(self):
        ...

    @abstractmethod
    def get_target_spec(self):
        ...

    @abstractmethod
    def rules(self):
        ...

    @property
    @abstractmethod
    def attempt(self):
        ...

    @property
    @abstractmethod
    def input(self):
        ...

    @property
    @abstractmethod
    def threads(self) -> int:
        ...

    @property
    @abstractmethod
    def resources(self):
        ...

    @property
    @abstractmethod
    def log(self):
        ...

    @abstractmethod
    def cleanup(self):
        ...

    @abstractmethod
    def get_wait_for_files(self):
        ...

    @abstractmethod
    def format_wildcards(self, string, **variables):
        ...

    @property
    @abstractmethod
    def needs_singularity(self):
        ...


class SingleJobExecutorInterface(ABC):
    @property
    @abstractmethod
    def rule(self):
        ...

    @abstractmethod
    def prepare(self):
        ...

    @property
    @abstractmethod
    def conda_env(self):
        ...

    @property
    @abstractmethod
    def container_img_path(self):
        ...

    @property
    @abstractmethod
    def env_modules(self):
        ...

    @property
    @abstractmethod
    def benchmark_repeats(self):
        ...

    @property
    @abstractmethod
    def benchmark(self):
        ...

    @property
    @abstractmethod
    def params(self):
        ...

    @property
    @abstractmethod
    def wildcards(self):
        ...

    @property
    @abstractmethod
    def shadow_dir(self):
        ...

    @property
    @abstractmethod
    def is_shadow(self):
        ...

    @property
    @abstractmethod
    def is_run(self):
        ...

    @property
    @abstractmethod
    def is_template_engine(self):
        ...

    @property
    @abstractmethod
    def message(self):
        ...


class GroupJobExecutorInterface(ABC):
    @property
    @abstractmethod
    def jobs(self):
        ...

    @property
    @abstractmethod
    def groupid(self):
        ...

    @property
    @abstractmethod
    def toposorted(self):
        ...


class DAGExecutorInterface(ABC):
    @abstractmethod
    def is_edit_notebook_job(self, job: ExecutorJobInterface):
        ...

    @abstractmethod
    def incomplete_external_jobid(self, job: ExecutorJobInterface):
        ...

    @abstractmethod
    def jobid(self, job: ExecutorJobInterface):
        ...

    @abstractmethod
    def get_sources(self):
        ...


class JobSchedulerExecutorInterface(ABC):
    @abstractmethod
    def executor_error_callback(self, exception):
        ...


class PersistenceExecutorInterface(ABC):
    @abstractmethod
    def cleanup(self):
        ...

    @property
    @abstractmethod
    def path(self):
        ...

    @property
    @abstractmethod
    def aux_path(self):
        ...


class WorkflowExecutorInterface(ABC):
    @property
    @abstractmethod
    def latency_wait(self) -> int:
        ...

    @property
    @abstractmethod
    def rerun_triggers(self) -> Optional[List[str]]:
        ...

    @property
    @abstractmethod
    def shadow_prefix(self) -> Optional[str]:
        ...

    @property
    @abstractmethod
    def conda_frontend(self) -> Optional[str]:
        ...

    @property
    @abstractmethod
    def conda_prefix(self) -> Optional[str]:
        ...

    @property
    @abstractmethod
    def conda_base_path(self) -> Optional[str]:
        ...

    @property
    @abstractmethod
    def singularity_args(self) -> Optional[str]:
        ...

    @property
    @abstractmethod
    def execute_subworkflows(self) -> bool:
        ...

    @property
    @abstractmethod
    def max_threads(self) -> Optional[int]:
        ...

    @property
    @abstractmethod
    def keep_metadata(self) -> bool:
        ...

    @property
    @abstractmethod
    def wrapper_prefix(self) -> Optional[str]:
        ...

    @property
    @abstractmethod
    def overwrite_threads(self) -> Dict[str, int]:
        ...

    @property
    @abstractmethod
    def overwrite_scatter(self) -> Dict[str, int]:
        ...

    @property
    @abstractmethod
    def local_groupid(self):
        ...

    @property
    @abstractmethod
    def conda_not_block_search_path_envvars(self):
        ...

    @property
    @abstractmethod
    def overwrite_configfiles(self):
        ...

    @property
    @abstractmethod
    def config_args(self):
        ...

    @property
    @abstractmethod
    def printshellcmds(self):
        ...

    @property
    @abstractmethod
    def scheduler_type(self):
        ...

    @property
    @abstractmethod
    def overwrite_resources(self):
        ...

    @property
    @abstractmethod
    def default_resources(self):
        ...

    @property
    @abstractmethod
    def overwrite_resource_scopes(self):
        ...

    @property
    @abstractmethod
    def resource_scopes(self):
        ...

    @abstractmethod
    def get_cache_mode(self, rule):
        ...

    @property
    @abstractmethod
    def output_file_cache(self):
        ...

    @property
    @abstractmethod
    def main_snakefile(self):
        ...

    @property
    @abstractmethod
    def persistence(self) -> PersistenceExecutorInterface:
        ...

    @property
    @abstractmethod
    def keep_metadata(self):
        ...

    @property
    @abstractmethod
    def linemaps(self):
        ...

    @property
    @abstractmethod
    def workdir_init(self):
        ...

    @property
    @abstractmethod
    def use_conda(self):
        ...

    @property
    @abstractmethod
    def use_singularity(self):
        ...

    @property
    @abstractmethod
    def use_env_modules(self):
        ...

    @property
    @abstractmethod
    def debug(self):
        ...

    @property
    @abstractmethod
    def cleanup_scripts(self):
        ...

    @property
    @abstractmethod
    def edit_notebook(self):
        ...

    @property
    @abstractmethod
    def sourcecache(self):
        ...

    @property
    @abstractmethod
    def verbose(self):
        ...

    @property
    @abstractmethod
    def jobscript(self):
        ...

    @property
    @abstractmethod
    def envvars(self):
        ...

    @property
    @abstractmethod
    def scheduler(self) -> JobSchedulerExecutorInterface:
        ...

    @property
    @abstractmethod
    def immediate_submit(self):
        ...

    @property
    @abstractmethod
    def default_remote_prefix(self):
        ...

    @property
    @abstractmethod
    def rules(self):
        ...

    @abstractmethod
    def get_rule(self, name):
        ...
