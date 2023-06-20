__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import importlib
import os
import pkgutil
from snakemake.logging import logger

executor_plugin_prefix = "snakemake_executor_"

# Executor plugins are externally installed plugins named "snakemake_executor_<name>"
# They should follow the same convention if on pip, snakemake-executor-<name>
executor_plugins = {
    name.replace(executor_plugin_prefix, ""): importlib.import_module(name)
    for _, name, _ in pkgutil.iter_modules()
    if name.startswith(executor_plugin_prefix)
}

# Valid Argument types (to distinguish from empty dataclasses)
ArgTypes = (str, int, float, bool)


def args_to_dataclass(args, dc):
    """
    Given snakemake args and a dataclass (dc) parse back into the Dataclass

    This helper function will select executor plugin namespaces arguments
    for a dataclass. It allows us to pass them from the custom executor ->
    custom argument parser -> back into dataclass -> snakemake.
    """
    # Iterate through the args, and parse those in the namespace
    prefix = get_plugin_name_from_dataclass(dc)
    kwargs = {}

    # Only derive those provided by the dataclass
    for arg in dir(dc):
        # Stick to those namespaced for the executor
        if not arg.startswith(prefix):
            continue
        value = getattr(args, arg)

        # And those that are in acceptable types
        if isinstance(value, ArgTypes):
            kwargs[arg] = value

    # At this point we want to convert the args <dataclasses> back into dataclass
    return dc(**kwargs)


def validate_executor_plugin(module):
    """
    Plugins must expose a common set of interfaces to be valid.
    """
    for name in ["add_args", "local_executor", "executor", "parse"]:
        if not hasattr(module, name):
            logger.warning(
                f"Executor plugin {module} is missing expected attribute {name}"
            )
            return False
    return True


def get_plugin_name_from_dataclass(dc):
    """
    Derive the module name (and executor name) from a dataclass.
    """
    return dc.__module__.split(".")[0].replace(executor_plugin_prefix, "")


def get_plugin_name(filepath):
    """
    Derive the name of a plugin from the filepath.

    This ensures that the plugin reported name (based on the module name)
    matches what we provide for the argument group.
    """
    return os.path.basename(
        os.path.dirname(filepath).replace(executor_plugin_prefix, "")
    )
