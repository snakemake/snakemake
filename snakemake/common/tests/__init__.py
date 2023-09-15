from abc import ABC, abstractmethod
from pathlib import Path
import shutil
from typing import List, Optional

import pytest
from snakemake import api, settings

from snakemake_interface_executor_plugins import ExecutorSettingsBase
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry


def handle_testcase(func):
    def wrapper(self, tmp_path):
        if self.expect_exception is None:
            try:
                return func(self, tmp_path)
            finally:
                self.cleanup_test()
        else:
            with pytest.raises(self.expect_exception):
                try:
                    return func(self, tmp_path)
                finally:
                    self.cleanup_test()

    return wrapper


class TestWorkflowsBase(ABC):
    __test__ = False
    expect_exception = None

    @abstractmethod
    def get_executor(self) -> str:
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

    def get_assume_shared_fs(self) -> bool:
        return True

    def get_envvars(self) -> List[str]:
        return []

    def cleanup_test(self):
        """This method is called after every testcase, also in case of exceptions.

        Override to clean up any test files (e.g. in remote storage).
        """
        pass

    def _run_workflow(self, test_name, tmp_path, deployment_method=frozenset()):
        test_path = Path(__file__).parent / "testcases" / test_name
        tmp_path = Path(tmp_path) / test_name
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
            storage_settings=settings.StorageSettings(
                default_remote_provider=self.get_default_remote_provider(),
                default_remote_prefix=self.get_default_remote_prefix(),
                assume_shared_fs=self.get_assume_shared_fs(),
            ),
            workdir=Path(tmp_path),
            snakefile=test_path / "Snakefile",
        )

        dag_api = workflow_api.dag(
            deployment_settings=settings.DeploymentSettings(
                deployment_method=deployment_method,
            ),
        )
        dag_api.execute_workflow(
            executor=self.get_executor(),
            executor_settings=self.get_executor_settings(),
            remote_execution_settings=settings.RemoteExecutionSettings(
                seconds_between_status_checks=0,
                envvars=self.get_envvars(),
            ),
        )
        snakemake_api.cleanup()

    @handle_testcase
    def test_simple_workflow(self, tmp_path):
        self._run_workflow("simple", tmp_path)

    @handle_testcase
    def test_group_workflow(self, tmp_path):
        self._run_workflow("groups", tmp_path)

    def _copy_test_files(self, test_path, tmp_path):
        shutil.copytree(test_path, tmp_path)

    def _common_settings(self):
        registry = ExecutorPluginRegistry()
        return registry.get_plugin(self.get_executor()).common_settings
