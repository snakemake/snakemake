from dataclasses import dataclass, fields
import hashlib
import os
import sys
from typing import Mapping, TypeVar, TYPE_CHECKING, Any
from snakemake_interface_executor_plugins.utils import format_cli_arg, join_cli_args
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry

from snakemake import common
from snakemake.settings import SharedFSUsage

if TYPE_CHECKING:
    from snakemake.workflow import Workflow

    TWorkflow = TypeVar("TWorkflow", bound="Workflow")
else:
    TWorkflow = Any


@dataclass
class SpawnedJobArgsFactory:
    workflow: TWorkflow

    def get_default_storage_provider_args(self) -> str:
        has_default_storage_provider = (
            self.workflow.storage_registry.default_storage_provider is not None
        )
        if has_default_storage_provider:
            return join_cli_args(
                [
                    format_cli_arg(
                        "--default-storage-prefix",
                        self.workflow.storage_settings.default_storage_prefix,
                    ),
                    format_cli_arg(
                        "--default-storage-provider",
                        self.workflow.storage_settings.default_storage_provider,
                    ),
                ]
            )
        else:
            return ""

    def _get_storage_provider_setting_items(self):
        for (
            plugin_name,
            tagged_settings,
        ) in self.workflow.storage_provider_settings.items():
            plugin = StoragePluginRegistry().get_plugin(plugin_name)
            if plugin.settings_cls:
                for field in fields(plugin.settings_cls):
                    unparse = field.metadata.get("unparse", lambda value: value)

                    def fmt_value(tag, value):
                        value = unparse(value)
                        if tag is not None:
                            return f"{tag}:{value}"
                        else:
                            return value

                    field_settings = [
                        fmt_value(tag, value)
                        for tag, value in tagged_settings.get_field_settings(
                            field.name
                        ).items()
                        if value is not None
                    ]
                    if field_settings:
                        yield plugin, field, field_settings

    def get_storage_provider_args(self):
        for plugin, field, field_settings in self._get_storage_provider_setting_items():
            if not field.metadata.get("env_var", False):
                cli_arg = plugin.get_cli_arg(field.name)
                yield format_cli_arg(cli_arg, field_settings)

    def get_storage_provider_envvars(self):
        return {
            plugin.get_envvar(field.name): " ".join(map(str, field_settings))
            for plugin, field, field_settings in self._get_storage_provider_setting_items()
            if "env_var" in field.metadata
        }

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
            "--set-resource-scopes",
            self.workflow.resource_settings.overwrite_resource_scopes,
        )

    def get_shared_fs_usage_arg(self):
        usage = self.workflow.storage_settings.shared_fs_usage
        return format_cli_arg(
            "--shared-fs-usage",
            usage if usage else "none",
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
            flag = f"--{query[-1].replace('_', '-')}"

        if invert and isinstance(value, bool):
            value = not value

        return format_cli_arg(flag, value, quote=quote)

    def envvars(self) -> Mapping[str, str]:
        envvars = {
            var: os.environ[var]
            for var in self.workflow.remote_execution_settings.envvars
        }
        envvars.update(self.get_storage_provider_envvars())
        return envvars

    def precommand(
        self,
        auto_deploy_default_storage_provider: bool = False,
        python_executable: str = "python",
    ) -> str:
        precommand = []
        if self.workflow.remote_execution_settings.precommand:
            precommand.append(self.workflow.remote_execution_settings.precommand)
        if (
            auto_deploy_default_storage_provider
            and self.workflow.storage_settings.default_storage_provider is not None
        ):
            package_name = StoragePluginRegistry().get_plugin_package_name(
                self.workflow.storage_settings.default_storage_provider
            )
            precommand.append(
                f"pip install --target '{common.PIP_DEPLOYMENTS_PATH}' {package_name}"
            )

        if (
            SharedFSUsage.SOURCES not in self.workflow.storage_settings.shared_fs_usage
            and self.workflow.remote_execution_settings.job_deploy_sources
        ):
            archive = self.workflow.source_archive
            default_storage_provider_args = self.get_default_storage_provider_args()
            storage_provider_args = " ".join(self.get_storage_provider_args())
            precommand.append(
                f"{python_executable} -m snakemake --deploy-sources "
                f"{archive.query} {archive.checksum} {default_storage_provider_args} "
                f"{storage_provider_args}"
            )

        return " && ".join(precommand)

    def general_args(
        self,
        pass_default_storage_provider_args: bool = True,
        pass_default_resources_args: bool = False,
    ) -> str:
        """Return a string to add to self.exec_job that includes additional
        arguments from the command line. This is currently used in the
        ClusterExecutor and CPUExecutor, as both were using the same
        code. Both have base class of the RealExecutor.
        """
        w2a = self.workflow_property_to_arg

        shared_deployment = (
            SharedFSUsage.SOFTWARE_DEPLOYMENT
            in self.workflow.storage_settings.shared_fs_usage
        )

        args = [
            "--force",
            "--target-files-omit-workdir-adjustment",
            "--keep-storage-local-copies",
            "--max-inventory-time 0",
            "--nocolor",
            "--notemp",
            "--no-hooks",
            "--nolock",
            "--ignore-incomplete",
            w2a("execution_settings.keep_incomplete"),
            w2a("output_settings.verbose"),
            w2a("rerun_triggers"),
            w2a(
                "execution_settings.cleanup_scripts",
                invert=True,
                flag="--skip-script-cleanup",
            ),
            w2a("execution_settings.shadow_prefix"),
            w2a("deployment_settings.deployment_method"),
            w2a("deployment_settings.conda_frontend"),
            w2a("deployment_settings.conda_prefix"),
            w2a(
                "conda_base_path",
                skip=not shared_deployment,
            ),
            w2a("deployment_settings.apptainer_prefix"),
            w2a("deployment_settings.apptainer_args"),
            w2a("resource_settings.max_threads"),
            self.get_shared_fs_usage_arg(),
            w2a(
                "execution_settings.keep_metadata", flag="--drop-metadata", invert=True
            ),
            w2a("workflow_settings.wrapper_prefix"),
            w2a("resource_settings.overwrite_threads", flag="--set-threads"),
            w2a("resource_settings.overwrite_scatter", flag="--set-scatter"),
            w2a("deployment_settings.conda_not_block_search_path_envvars"),
            w2a("overwrite_configfiles", flag="--configfiles"),
            w2a("config_settings.config_args", flag="--config"),
            w2a("output_settings.printshellcmds"),
            w2a("execution_settings.latency_wait"),
            w2a("scheduling_settings.scheduler", flag="--scheduler"),
            format_cli_arg(
                "--scheduler-solver-path",
                os.path.dirname(sys.executable),
                skip=not shared_deployment,
            ),
            w2a(
                "overwrite_workdir",
                flag="--directory",
                skip=self.workflow.storage_settings.assume_common_workdir,
            ),
            self.get_set_resources_args(),
            self.get_resource_scopes_args(),
        ]
        args.extend(self.get_storage_provider_args())
        if pass_default_storage_provider_args:
            args.append(self.get_default_storage_provider_args())
        if pass_default_resources_args:
            args.append(w2a("resource_settings.default_resources", attr="args"))

        return join_cli_args(args)
