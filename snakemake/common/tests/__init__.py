from abc import ABC, abstractmethod
from pathlib import Path
import shutil
from typing import Optional
from snakemake import api, settings

from snakemake_interface_executor_plugins import ExecutorSettingsBase
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry

class TestWorkflow(ABC):

    @abstractmethod
    def get_executor(self):
        ...

    @abstractmethod
    def get_executor_settings(self) -> Optional[ExecutorSettingsBase]:
        ...

    @abstractmethod
    def get_default_remote_provider(self) -> Optional[str]:
        ...
    
    @abstractmethod
    def get_default_remote_prefix(self) -> Optional[str]:
        ...

    def _run_workflow(self, test_name, tmp_path, deployment_method=frozenset()):
        test_path = Path(__file__).parent / "testcases" / test_name
        self._copy_test_files(test_path, tmp_path)

        if self._common_settings().local_exec:
            cores = 3
            nodes = None
        else:
            cores = 1
            nodes = 3

        snakemake_api = api.SnakemakeApi(
            settings.OutputSettings(
                verbose=True,
            ),
        )
        workflow_api = snakemake_api.workflow(
            resource_settings=settings.ResourceSettings(
                cores=cores,
                nodes=nodes,
            ),
            workdir=Path(tmp_path),
        )

        dag_api = workflow_api.dag(
            deployment_settings=settings.DeploymentSettings(
                deployment_method=deployment_method,
            ),
        )
        dag_api.execute_workflow(
            executor=self.get_executor(),
            executor_settings=self.get_executor_settings(),
        )
        snakemake_api.cleanup()

    def test_simple_workflow(self, tmp_path):
        self._run_workflow("simple_workflow", tmp_path)
    
    def test_group_workflow(self, tmp_path):
        self._run_workflow("group_workflow", tmp_path)

    def _copy_test_files(self, test_path, tmp_path):
        shutil.copytree(test_path, tmp_path)
    
    def _common_settings(self):
        registry = ExecutorPluginRegistry()
        return registry[self.get_executor()].common_settings