from dataclasses import dataclass, fields
import hashlib
from itertools import chain
import os
import sys
from typing import Callable, Mapping, TypeVar, TYPE_CHECKING, Any
from snakemake_interface_executor_plugins.utils import format_cli_arg, join_cli_args
from snakemake_interface_executor_plugins.settings import CommonSettings
from snakemake.resources import ParsedResource
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry

from snakemake import PIP_DEPLOYMENTS_PATH
from snakemake.io import get_flag_value, is_flagged
from snakemake.settings.types import SharedFSUsage
from snakemake_interface_common.exceptions import WorkflowError

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
                        if callable(value):
                            raise WorkflowError(
                                f"Invalid setting for plugin {plugin_name}. Unable "
                                "to pass callable value as a setting to spawned jobs."
                            )
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
                if isinstance(field_settings[0], bool) and field.default is not True:
                    # so far no tagged settings for flags with default value False
                    assert len(field_settings) == 1
                    field_settings = field_settings[0]
                yield format_cli_arg(cli_arg, field_settings)

    def get_storage_provider_envvars(self):
        return {
            plugin.get_envvar(field.name): " ".join(map(str, field_settings))
            for plugin, field, field_settings in self._get_storage_provider_setting_items()
            if field.metadata.get("env_var")
        }

    def get_set_resources_args(self):
        def get_orig_arg(value):
            if isinstance(value, ParsedResource):
                return value.orig_arg
            else:
                return value

        return [
            format_cli_arg(
                "--set-resources",
                [
                    f"{rule}:{name}={get_orig_arg(value)}"
                    for rule, res in self.workflow.resource_settings.overwrite_resources.items()
                    for name, value in res.items()
                ],
                skip=not self.workflow.resource_settings.overwrite_resources,
                base64_encode=True,
            ),
            format_cli_arg(
                "--set-threads",
                [
                    f"{rule}={get_orig_arg(value)}"
                    for rule, value in self.workflow.resource_settings.overwrite_threads.items()
                ],
                skip=not self.workflow.resource_settings.overwrite_threads,
                base64_encode=True,
            ),
        ]

    def get_resource_scopes_arg(self):
        return format_cli_arg(
            "--set-resource-scopes",
            self.workflow.resource_settings.overwrite_resource_scopes,
        )

    def get_shared_fs_usage_arg(self, executor_common_settings: CommonSettings):
        if executor_common_settings.spawned_jobs_assume_shared_fs:
            usage = SharedFSUsage.all()
        else:
            usage = self.workflow.storage_settings.shared_fs_usage
        return format_cli_arg(
            "--shared-fs-usage",
            usage if usage else "none",
        )

    def get_group_args(self):
        group_settings = self.workflow.group_settings
        groups = format_cli_arg(
            "--groups",
            [
                f"{rule}={group}"
                for rule, group in group_settings.overwrite_groups.items()
            ],
            skip=not group_settings.overwrite_groups,
        )
        group_components = format_cli_arg(
            "--group-components",
            [f"{group}={n}" for group, n in group_settings.group_components.items()],
            skip=not group_settings.group_components,
        )
        return join_cli_args([groups, group_components])

    def workflow_property_to_arg(
        self,
        property,
        flag=None,
        base64_encode=False,
        skip=False,
        force=False,
        invert=False,
        attr=None,
        convert_value: Callable = None,
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

        if force:
            assert isinstance(value, bool)
            value = True

        if convert_value is not None and value is not None:
            value = convert_value(value)

        return format_cli_arg(flag, value, base64_encode=base64_encode)

    def envvars(self) -> Mapping[str, str]:
        assert self.workflow.remote_execution_settings is not None
        envvars = {
            var: os.environ[var]
            for var in chain(
                self.workflow.remote_execution_settings.envvars, self.workflow.envvars
            )
        }
        envvars.update(self.get_storage_provider_envvars())
        return envvars

    def precommand(
        self,
        executor_common_settings: CommonSettings,
        python_executable: str = "python",
    ) -> str:
        precommand = []
        if self.workflow.remote_execution_settings.precommand:
            precommand.append(self.workflow.remote_execution_settings.precommand)
        if (
            executor_common_settings.auto_deploy_default_storage_provider
            and self.workflow.storage_settings.default_storage_provider is not None
        ):
            packages_to_install = set(
                StoragePluginRegistry().get_plugin_package_name(pkg)
                for pkg in self.workflow.storage_provider_settings.keys()
            )
            pkgs = " ".join(packages_to_install)

            precommand.append(f"pip install --target '{PIP_DEPLOYMENTS_PATH}' {pkgs}")

        if (
            SharedFSUsage.SOURCES not in self.workflow.storage_settings.shared_fs_usage
            and self.workflow.remote_execution_settings.job_deploy_sources
            and not executor_common_settings.can_transfer_local_files
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

    def get_configfiles_arg(self):
        # If not shared FS, use relpath for the configfiles (deployed via source archive)
        # Source archive creation ensures that configfiles are in a subdir of the CWD
        # and errors otherwise.
        if SharedFSUsage.SOURCES not in self.workflow.storage_settings.shared_fs_usage:
            configfiles = [
                os.path.relpath(f) for f in self.workflow.overwrite_configfiles
            ]
        else:
            configfiles = self.workflow.overwrite_configfiles
        if configfiles:
            return format_cli_arg("--configfiles", configfiles)
        else:
            return ""

    def general_args(
        self,
        executor_common_settings: CommonSettings,
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

        # base64 encode the prefix to ensure that eventually unexpanded env vars
        # are not replaced with values (or become empty if missing) by the shell
        if executor_common_settings.non_local_exec and not self.workflow.remote_exec:
            # this is used when the main process submits a job via a remote executor
            local_storage_prefix = w2a(
                "storage_settings.remote_job_local_storage_prefix",
                flag="--local-storage-prefix",
                base64_encode=True,
            )
        else:
            # this is used in the main process when submitting jobs to a local executor
            # or when a remote executor spawns an inner executor
            local_storage_prefix = w2a(
                "storage_settings.local_storage_prefix", base64_encode=True
            )

        args = [
            "--force",
            "--target-files-omit-workdir-adjustment",
            w2a(
                "storage_settings.keep_storage_local",
                flag="--keep-storage-local-copies",
            ),
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
                "storage_settings.wait_for_free_local_storage",
                convert_value="{}s".format,
            ),
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
            w2a("deployment_settings.apptainer_args", base64_encode=True),
            w2a("resource_settings.max_threads"),
            self.get_shared_fs_usage_arg(executor_common_settings),
            w2a(
                "execution_settings.keep_metadata", flag="--drop-metadata", invert=True
            ),
            w2a("workflow_settings.wrapper_prefix"),
            w2a("resource_settings.overwrite_scatter", flag="--set-scatter"),
            w2a("deployment_settings.conda_not_block_search_path_envvars"),
            w2a("config_settings.config_args", flag="--config"),
            w2a("output_settings.printshellcmds"),
            w2a("output_settings.benchmark_extended"),
            w2a("execution_settings.latency_wait"),
            w2a("scheduling_settings.scheduler", flag="--scheduler"),
            w2a("workflow_settings.cache"),
            local_storage_prefix,
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
            self.get_resource_scopes_arg(),
            self.get_configfiles_arg(),
        ]
        args.extend(self.get_storage_provider_args())
        args.extend(self.get_set_resources_args())
        if executor_common_settings.pass_default_storage_provider_args:
            args.append(self.get_default_storage_provider_args())
        if executor_common_settings.pass_default_resources_args:
            args.append(
                w2a(
                    "resource_settings.default_resources",
                    attr="args",
                    base64_encode=True,
                )
            )
        if executor_common_settings.pass_group_args:
            args.append(self.get_group_args())

        from snakemake.logging import logger

        logger.debug(f"General args: {args}")

        return join_cli_args(args)
