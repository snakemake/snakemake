import importlib
from dataclasses import fields, MISSING
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List

import click

from snakemake_interface_common._common import dataclass_field_to_argument_args
from snakemake_interface_common.plugin_registry.plugin import PluginBase
from snakemake_interface_common.plugin_registry import PluginRegistryBase
from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
from snakemake_interface_storage_plugins.registry import StoragePluginRegistry
from snakemake_interface_report_plugins.registry import ReportPluginRegistry
from snakemake_interface_logger_plugins.registry import LoggerPluginRegistry
from snakemake_interface_scheduler_plugins.registry import SchedulerPluginRegistry

PLUGIN_REGISTRIES: Dict[str, type[PluginRegistryBase]] = {
    "executor": ExecutorPluginRegistry,
    "storage": StoragePluginRegistry,
    "report": ReportPluginRegistry,
    "logger": LoggerPluginRegistry,
    "scheduler": SchedulerPluginRegistry,
}


def get_registry(plugin_type: str) -> PluginRegistryBase:
    """Instantiate a plugin registry by type name.

    Args:
        plugin_type: One of "executor", "storage", "report", "logger", "scheduler".

    Raises:
        KeyError: If plugin_type is not a recognized type.
    """
    return PLUGIN_REGISTRIES[plugin_type]()


def get_install_path(plugin_name: str, module_prefix: str) -> str:
    """Get the filesystem path where a plugin package is installed.

    Returns ``'built-in'`` for plugins bundled with snakemake, or
    ``'unknown'`` if the path cannot be determined.
    """
    package_name = f"{module_prefix}{plugin_name.replace('-', '_')}"
    try:
        mod = importlib.import_module(package_name)
        path = getattr(mod, "__file__", None)
        if path:
            return str(Path(path).parent)
    except ImportError:
        return "built-in"

    return "unknown"


def click_params_from_plugin(plugin: PluginBase) -> List[click.Parameter]:
    """Convert a plugin's settings dataclass into Click options.

    Uses ``snakemake_interface_common._common.dataclass_field_to_argument_args``
    to build the ``(args, kwargs)`` pairs for each field, then applies the same
    post-processing that
    ``snakemake_interface_common.plugin_registry.plugin.PluginBase.register_cli_args``
    does (env_var, metavar, tagged values, parse_func, Path help text) before
    converting to ``click.Option`` instances.
    """
    if not plugin.has_settings_cls():
        return []

    params: List[click.Parameter] = []
    for thefield in fields(plugin.settings_cls):  # type: ignore
        prefixed_name = plugin._get_prefixed_name(thefield.name).replace("-", "_")
        args, kwargs = dataclass_field_to_argument_args(thefield, name=prefixed_name)

        if thefield.metadata.get("env_var"):
            kwargs["env_var"] = plugin.get_envvar(thefield.name)

        if "metavar" not in kwargs:
            kwargs["metavar"] = "VALUE"

        if thefield.type == Path:
            kwargs["help"] = kwargs.get("help", "") + (
                " User dir (~) and environment variables are properly interpreted."
            )

        if thefield.type == bool:
            kwargs.pop("metavar", None)
        elif plugin.support_tagged_values:
            if thefield.metadata.get("nargs", None) is not None:
                raise ValueError(
                    f"Plugin {plugin.name} supports tagged values but specifies args "
                    "with multiple values in its settings class."
                )
            kwargs["nargs"] = "+"
            kwargs["type"] = str
            kwargs["help"] = kwargs.get("help", "") + (
                " Can be specified multiple times to set different "
                "values for different tags."
            )
            kwargs["metavar"] = f"[TAG::]{kwargs['metavar']}"

        if thefield.metadata.get("parse_func"):
            kwargs["type"] = str

        name = args[-1]  # e.g. --executor-slurm-partition
        click_kwargs: Dict[str, Any] = {"help": kwargs.get("help", "")}

        if thefield.type == bool:
            click_kwargs["is_flag"] = True
            default = kwargs.get("default", MISSING)
            click_kwargs["default"] = default if default is not MISSING else False
        else:
            default = kwargs.get("default", MISSING)
            click_kwargs["default"] = default if default is not MISSING else None

            if "choices" in kwargs:
                click_kwargs["type"] = click.Choice(kwargs["choices"])
            elif kwargs.get("type") not in (None, str):
                click_kwargs["type"] = kwargs["type"]

            if kwargs.get("nargs") == "+":
                click_kwargs["multiple"] = True

            if kwargs.get("required"):
                click_kwargs["required"] = True

        if "env_var" in kwargs:
            click_kwargs["envvar"] = kwargs["env_var"]

        if "metavar" in kwargs and not click_kwargs.get("is_flag"):
            click_kwargs["metavar"] = kwargs["metavar"]

        params.append(click.Option([name], **click_kwargs))

    return params


def click_params_from_registry(
    registry: PluginRegistryBase,
) -> List[click.Parameter]:
    """Get Click params for all plugins in a registry."""
    params: List[click.Parameter] = []
    for _, plugin in registry.plugins.items():
        params.extend(click_params_from_plugin(plugin))
    return params


def add_plugin_options(*plugin_types: str):
    """Decorator that attaches plugin options to a Click command.

    Args:
        plugin_types: Plugin type names (e.g. ``"executor"``, ``"storage"``).

    Usage::

        @click.command()
        @workflow_options
        @add_plugin_options("executor", "storage")
        def run(snakefile, directory, **kwargs):
            ...
    """

    def decorator(func):
        for plugin_type in plugin_types:
            registry = get_registry(plugin_type)
            for param in click_params_from_registry(registry):
                func.params.append(param)
        return func

    return decorator


def plugin_settings_from_kwargs(plugin: PluginBase, kwargs: Dict[str, Any]) -> Any:
    """Extract plugin settings from a flat Click kwargs dict.

    Builds a ``SimpleNamespace`` that mimics an argparse ``Namespace``
    so that ``PluginBase.get_settings()`` works unchanged.
    """
    ns = SimpleNamespace(**kwargs)
    return plugin.get_settings(ns)
