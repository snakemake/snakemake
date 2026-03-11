"""Shared context for subcommands.

SnakemakeContext accumulates settings and provides managed context managers
that handle SnakemakeApi lifecycle. Subcommands construct a context from
their click options, then call .workflow() or .dag() to get the API objects
they need.

"""

from contextlib import contextmanager
from pathlib import Path
from typing import Optional

from snakemake.api import SnakemakeApi
from snakemake.settings.types import (
    ConfigSettings,
    DAGSettings,
    DeploymentSettings,
    OutputSettings,
    ResourceSettings,
    StorageSettings,
    WorkflowSettings,
)


class SnakemakeContext:
    """Bridges click subcommands and the Snakemake API.

    Holds all settings needed to construct API objects. Context managers
    handle SnakemakeApi enter/exit so callers get proper cleanup.
    """

    def __init__(
        self,
        snakefile: Optional[str] = None,
        workdir: Optional[str] = None,
        output_settings: Optional[OutputSettings] = None,
        resource_settings: Optional[ResourceSettings] = None,
        config_settings: Optional[ConfigSettings] = None,
        storage_settings: Optional[StorageSettings] = None,
        workflow_settings: Optional[WorkflowSettings] = None,
        deployment_settings: Optional[DeploymentSettings] = None,
        dag_settings: Optional[DAGSettings] = None,
    ):
        self.snakefile = Path(snakefile) if snakefile else None
        self.workdir = Path(workdir) if workdir else None
        self.output_settings = output_settings or OutputSettings()
        self.resource_settings = resource_settings or ResourceSettings()
        self.config_settings = config_settings or ConfigSettings()
        self.storage_settings = storage_settings or StorageSettings()
        self.workflow_settings = workflow_settings or WorkflowSettings()
        self.deployment_settings = deployment_settings or DeploymentSettings()
        self.dag_settings = dag_settings or DAGSettings()

    @contextmanager
    def workflow(self):
        """Yield a WorkflowApi inside a managed SnakemakeApi context."""
        with SnakemakeApi(self.output_settings) as api:
            workflow_api = api.workflow(
                resource_settings=self.resource_settings,
                config_settings=self.config_settings,
                storage_settings=self.storage_settings,
                workflow_settings=self.workflow_settings,
                deployment_settings=self.deployment_settings,
                snakefile=self.snakefile,
                workdir=self.workdir,
            )
            yield workflow_api

    @contextmanager
    def dag(self):
        """Yield a DAGApi inside a managed SnakemakeApi context."""
        with self.workflow() as workflow_api:
            yield workflow_api.dag(dag_settings=self.dag_settings)
