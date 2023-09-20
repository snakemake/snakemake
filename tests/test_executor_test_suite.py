from pathlib import Path
from typing import Optional
import snakemake.common.tests

from snakemake_interface_executor_plugins import ExecutorSettingsBase
from snakemake_executor_plugin_cluster_generic import ExecutorSettings


class TestWorkflows(snakemake.common.tests.TestWorkflowsBase):
    __test__ = True

    def get_executor(self) -> str:
        return "cluster-generic"

    def _get_cmd(self) -> str:
        return str((Path(__file__).parent / "test_group_jobs" / "qsub").absolute())

    def get_executor_settings(self) -> Optional[ExecutorSettingsBase]:
        return ExecutorSettings(submit_cmd=self._get_cmd())

    def get_default_remote_provider(self) -> Optional[str]:
        return None

    def get_default_remote_prefix(self) -> Optional[str]:
        return None
