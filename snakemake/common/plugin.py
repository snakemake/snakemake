__author__ = "Hielke Walinga"
__copyright__ = "Copyright 2022, Hielke Walinga"
__email__ = "hielkewalinga@gmail.com"
__license__ = "MIT"

from os import path
import sys
import logger
import importlib
import pkgutil
from functools import lru_cache
import collections

from snakemake.exceptions import PluginException


def find_plugins(prefix):
    """Iterator to find all plugins by their prefix of the package name."""
    for finder, name, ispkg in pkgutil.iter_modules():
        # Check if the *package name* of the module
        # start with the correct prefix.
        if cached_packages_distributions()[name][0].startswith(prefix):
            yield importlib.import_module(name)


plugin_required_mod_attrs = [
    "__version__", "__author__", "__copyright__", "__email__", "__license__"
]


def verify_plugin(
    mod,
    mod_attrs=plugin_required_mod_attrs,
    extra_attrs=[],
    checks={}
):
    """Confirms a module (just an object) has the required attributes

    Raises a PluginException if this is not true.

    Also takes a checks dictionary with a function and an optional exception.
    The check should return None or False to pass.
    """
    package_name = cached_packages_distributions()[mod.__name__][0]
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


def load_plugins(plugin_modules, globals_dict, extra_attrs, checks):
    """Load plugins into the globals dict, yielding all the modules in the process.
    Additionally some checks are performed on the plugins.
    NB. Don't forget to exhaust the created generated to have an effect."""
    for plugin_module in plugin_modules:
        plugin_name = plugin_module.__name__
        package_name = cached_packages_distributions()[plugin_name][0]
        if plugin_module.__name__ in globals_dict:
            logger.debug("Plugin {package_name} already loaded.")
            continue

        try:
            verify_plugin(plugin_module, extra_attrs=extra_attrs, checks=checks)
        except PluginException as e:
            logger.warning(f"Plugin {package_name} incorrect: {e}")
            continue

        yield plugin_module
        snakemake_submodule_name = (
            getattr(plugin_module, "snakemake_submodule_name") or plugin_name
        )
        globals_dict[snakemake_submodule_name] = plugin_module
        globals_dict[snakemake_submodule_name.lower()] = plugin_module


@lru_cache(max_size=None)
def cached_packages_distributions():
    """Return a mapping of top-level packages to their distributions."""

    # Source: https://github.com/python/cpython/blob/9006b4471cc4d74d0cbf782d841a330b0f46a3d0/Lib/importlib/metadata/__init__.py#L1031

    # TODO: Transition instructions:
    # importlib.metadata becomes part of stdlib from 3.8 onwards
    # importlib.metadata will have package_distributions from 3.10 onwards

    try:
        from importlib.metadata import packages_distributions as pkgs_dist
        return pkgs_dist()
    except ImportError:
        ...  # We are older than 3.10, import distributions instead.

    if hasattr(importlib, 'metadata'):
        from importlib.metadata import distributions
    else:
        # We are older than 3.8,
        # import from external importlib_metadata instead.
        from importlib_metadata import distributions

    pkg_to_dist = collections.defaultdict(list)
    for dist in distributions():
        for pkg in (dist.read_text('top_level.txt') or '').split():
            pkg_to_dist[pkg].append(dist.metadata['Name'])
    return dict(pkg_to_dist)


def internal_submodules(paths):
    for remote_submodule in pkgutil.iter_modules(paths):
        module_name = remote_submodule.name

        if module_name == 'common':  # Does not have a specific RemoteProvider
            continue

        module_path = (
            path.join(remote_submodule.module_finder.path, module_name)
            + ".py"
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
