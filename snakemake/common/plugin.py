__author__ = "Hielke Walinga"
__copyright__ = "Copyright 2022, Hielke Walinga"
__email__ = "hielkewalinga@gmail.com"
__license__ = "MIT"

import importlib
import pkgutil
import sys
from os import path

from snakemake.common.future import packages_distributions
from snakemake.exceptions import PluginException
from snakemake.logging import logger

cached_packages_distributions = packages_distributions()


def find_plugins(prefix):
    """Iterator to find all plugins by their prefix of the package name."""
    for finder, name, ispkg in pkgutil.iter_modules():
        # Check if the *package name* of the module
        # start with the correct prefix.
        try:
            print("--")
            print(name)
            packages = cached_packages_distributions[name]
            print(packages)
        except KeyError:
            continue
        if packages[0].startswith(prefix):
            yield importlib.import_module(name)


plugin_required_mod_attrs = [
    "__version__",
    "__author__",
    "__copyright__",
    "__email__",
    "__license__",
    "__min_snakemake_version__",
]


def verify_plugin(
    mod, mod_attrs=plugin_required_mod_attrs, extra_attrs=[], checks={}, prefix=None
):
    """Confirms a module (just an object) has the required attributes

    Raises a PluginException if this is not true.

    Also takes a checks dictionary with a function and an optional exception.
    The check should return None or False to pass.
    """
    package_name = cached_packages_distributions[mod.__name__][0]
    if prefix:
        if not package_name.startswith(prefix):
            raise PluginException(
                f"Plugin {package_name} is not distributed "
                f"with a package starting with {prefix}."
            )
    for a in mod_attrs + extra_attrs:
        if not hasattr(mod, a):
            raise PluginException(
                f"Plugin {package_name} does not provide attribute {a}."
            )
    for check, exception in checks.items():
        if check(mod):
            raise (
                exception
                or PluginException(
                    f"Plugin {package_name} fails check {check.__name__}."
                )
            )


def load_plugins(plugin_modules, globals_dict, verify_func, alternative_mod_name=None):
    """Load plugins into the globals dict, yielding all the modules in the process.
    Additionally some checks are performed on the plugins.
    NB. Don't forget to exhaust the created generated to have an effect."""
    for plugin_module in plugin_modules:
        plugin_name = plugin_module.__name__
        package_name = cached_packages_distributions[plugin_name][0]
        if plugin_module.__name__ in globals_dict:
            logger.debug("Plugin {package_name} already loaded.")
            continue

        try:
            verify_func(plugin_module)
        except PluginException as e:
            logger.warning(f"Plugin {package_name} incorrect: {e}")
            continue

        yield plugin_module
        snakemake_submodule_name = getattr(
            plugin_module,
            alternative_mod_name or "snakemake_submodule_name",
            plugin_name,
        )
        globals_dict[snakemake_submodule_name] = plugin_module
        globals_dict[snakemake_submodule_name.lower()] = plugin_module


def internal_submodules(paths):
    for remote_submodule in pkgutil.iter_modules(paths):
        module_name = remote_submodule.name

        if module_name == "common":  # Does not have a specific RemoteProvider
            continue

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

        yield module
