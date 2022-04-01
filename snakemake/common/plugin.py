__author__ = "Hielke Walinga"
__copyright__ = "Copyright 2022, Hielke Walinga"
__email__ = "hielkewalinga@gmail.com"
__license__ = "MIT"

import importlib
import pkgutil
import sys
from os import path

# TODO: When 3.9 deprecated, import from importlib.metadata
from snakemake.common.future import packages_distributions
from snakemake.exceptions import PluginException
from snakemake.logging import logger
from snakemake.utils import min_version

cached_packages_distributions = packages_distributions()


def find_plugins(prefix):
    """Iterator to find all plugins by their prefix of the package name."""
    plugins = []
    for finder, name, ispkg in pkgutil.iter_modules():
        # Check if the *package name* of the module
        # start with the correct prefix.
        try:
            packages = cached_packages_distributions[name]
        except KeyError:
            continue
        for p in packages[1:]:
            if p.startswith(prefix):
                raise PluginException(
                    "A plugin is installed onto a module namespace taken. "
                    f"These packages conflict: {packages}."
                )
        if packages[0].startswith(prefix):
            plugins.append(importlib.import_module(name))

    return plugins


plugin_required_mod_attrs = [
    "__version__",
    "__author__",
    "__copyright__",
    "__email__",
    "__license__",
    "__min_snakemake_version__",
]


def verify_plugin(
    mod, mod_attrs=plugin_required_mod_attrs, extra_attrs=[], prefix=None
):
    """Confirms a module (just an object) has the required attributes

    Raises a PluginException if this is not true.

    Also takes a checks dictionary with a function and an optional exception.
    The check should return None or False to pass.
    """
    package_name = cached_packages_distributions[mod.__name__][0]
    if prefix and not package_name.startswith(prefix):
        raise PluginException(
            f"Plugin {package_name} is not distributed "
            f"with a package starting with {prefix}."
        )
    for a in mod_attrs + extra_attrs:
        if not hasattr(mod, a):
            raise PluginException(
                f"Plugin {package_name} does not provide attribute {a}."
            )
    try:
        min_version(mod.__min_snakemake_version__)
    except Exception:
        raise PluginException(
            "Plugin {mod.__name__} needs minimum Snakemake version {mod.__min_snakemake_version__}."
        )


def load_plugin(plugin_module, globals_dict, verify_func, alternative_mod_name=None):
    """Load plugins into the globals dict.
    Additionally some checks are performed on the plugins.
    The module is bround with the name and the lowercase of the name.

    plugin_module: module of the plugin
    globals_dict: result of globals to the set module on
    verify_func: function that verifies plugin or throw PluginException if wrong
    alternative_mod_name: source name on the plugin module that provides the name
        to assign the module on the globals_dict.

    returns True if plugin is loaded.
    """
    plugin_name = plugin_module.__name__
    package_name = cached_packages_distributions[plugin_name][0]
    if plugin_module.__name__ in globals_dict:
        logger.debug("Plugin {package_name} already loaded.")
        return True

    try:
        verify_func(plugin_module)
    except PluginException as e:
        logger.warning(f"Plugin {package_name} incorrect: {e}")
        return False

    snakemake_submodule_name = getattr(
        plugin_module,
        alternative_mod_name or "snakemake_submodule_name",
        plugin_name,
    )
    globals_dict[snakemake_submodule_name] = plugin_module
    globals_dict[snakemake_submodule_name.lower()] = plugin_module
    return True


def internal_submodules(paths):
    """Loads all modules found on the paths and returns them as a list."""
    submodules = []
    for remote_submodule in pkgutil.iter_modules(paths):
        module_name = remote_submodule.name

        module_path = (
            path.join(remote_submodule.module_finder.path, module_name) + ".py"
        )

        spec = importlib.util.spec_from_file_location(module_name, module_path)
        module = importlib.util.module_from_spec(spec)

        try:
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
        except Exception as e:
            logger.debug(f"Autoloading {module_name} failed: {e}")
            continue

        submodules.append(module)

    return submodules
