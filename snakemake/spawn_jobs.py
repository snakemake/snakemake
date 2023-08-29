from dataclasses import dataclass
import os
import sys
from snakemake_interface_executor_plugins.utils import format_cli_arg, join_cli_args

from snakemake.workflow import Workflow

@dataclass
class SpawnedJobArgsFactory:
    workflow: Workflow

    def get_default_remote_provider_args(self):
        return join_cli_args(
            [
                format_cli_arg("--default-remote-prefix", self.workflow.storage_settings.default_remote_prefix),
                format_cli_arg("--default-remote-provider", self.workflow.storage_settings.default_remote_provider.name),
            ]
        )

    def get_set_resources_args(self):
        return format_cli_arg(
            "--set-resources",
            [
                f"{rule}:{name}={value}"
                for rule, res in self.workflow.resource_settings.overwrite_resources.items()
                for name, value in res.items()
            ],
            skip=not self.workflow.resource_settings.overwrite_resources,
        )

    def get_resource_scopes_args(self):
        return format_cli_arg(
            "--set-resource-scopes", self.workflow.overwrite_resource_scopes
        )

    def workflow_property_to_arg(
        self, property, flag=None, quote=True, skip=False, invert=False, attr=None
    ):
        if skip:
            return ""

        # Get the value of the property. If property is nested, follow the hierarchy until
        # reaching the final value.
        query = property.split(".")
        base = self.workflow
        for prop in query[:-1]:
            base = getattr(base, prop)
        value = getattr(base, query[-1])

        if value is not None and attr is not None:
            value = getattr(value, attr)

        if flag is None:
            flag = f"--{property.replace('_', '-')}"

        if invert and isinstance(value, bool):
            value = not value

        return format_cli_arg(flag, value, quote=quote)

    def general_args(
            self,
            pass_default_remote_provider_args: bool=True, 
            pass_default_resources_args: bool=False
        ):
        """Return a string to add to self.exec_job that includes additional
        arguments from the command line. This is currently used in the
        ClusterExecutor and CPUExecutor, as both were using the same
        code. Both have base class of the RealExecutor.
        """
        w2a = self.workflow_property_to_arg

        args = [
                "--force",
                "--keep-target-files",
                "--keep-remote",
                "--max-inventory-time 0",
                "--nocolor",
                "--notemp",
                "--no-hooks",
                "--nolock",
                "--ignore-incomplete",
                w2a("execution_settings.keep_incomplete"),
                w2a("rerun_triggers"),
                w2a("execution_settings.cleanup_scripts", flag="--skip-script-cleanup"),
                w2a("execution_settings.shadow_prefix"),
                w2a("deployment_settings.deployment_method"),
                w2a("deployment_settings.conda_frontend"),
                w2a("deployment_settings.conda_prefix"),
                w2a("conda_base_path", skip=not self.assume_shared_fs),
                w2a("deployment_settings.use_singularity"),
                w2a("deployment_settings.apptainer_prefix"),
                w2a("deployment_settings.apptainer_args"),
                w2a("execution_settings.execute_subworkflows", flag="--no-subworkflows", invert=True),
                w2a("resource_settings.max_threads"),
                w2a("deployment_settings.use_env_modules", flag="--use-envmodules"),
                w2a("execution_settings.keep_metadata", flag="--drop-metadata", invert=True),
                w2a("execution_settings.wrapper_prefix"),
                w2a("resource_settings.overwrite_threads", flag="--set-threads"),
                w2a("overwrite_scatter", flag="--set-scatter"),
                w2a("deployment_settings.conda_not_block_search_path_envvars"),
                w2a("overwrite_configfiles", flag="--configfiles"),
                w2a("config_settings.config_args", flag="--config"),
                w2a("output_settings.printshellcmds"),
                w2a("latency_wait"),
                w2a("scheduler_settings.scheduler_type", flag="--scheduler"),
                format_cli_arg(
                    "--scheduler-solver-path",
                    os.path.dirname(sys.executable),
                    skip=not self.workflow.storage_settings.assume_shared_fs,
                ),
                self.get_set_resources_args(),
                self.get_resource_scopes_args(),
            ]
        if pass_default_remote_provider_args:
            args.append(self.get_default_remote_provider_args())
        if pass_default_resources_args:
            args.append(w2a("resource_settings.default_resources", attr="args"))

        return join_cli_args(args)