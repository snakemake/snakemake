__author__ = "Hielke Walinga"
__copyright__ = "Copyright 2022, Hielke Walinga"
__email__ = "hielkewalinga@gmail.com"
__license__ = "MIT"

import importlib
import pkgutil
from functools import lru_cache
import collections

from snakemake.exceptions import PluginException


def find_plugins(prefix):
    """Iterator to find all plugins by their prefix of the package name."""
    for finder, name, ispkg in pkgutil.iter_modules():
        # Check if the *package name* of the module start with the correct prefix.
        if cached_packages_distributions()[name][0].startswith(prefix):
            yield importlib.import_module(name)


def confirm_mod_has_correct_attrs(
    mod,
    mod_attrs=["__version__", "__copyright__", "__email__", "__licence__"],
    extra_attrs=[],
):
    """Confirms a module (just an object) has the required attributes

    Raises a PluginException if this is not true.
    """
    for a in mod_attrs + extra_attrs:
        if not hasattr(mod, a):
            package_name = cached_packages_distributions()[mod.__name__][0]
            raise PluginException(f"Plugin {package_name} does not provide attribute {a}.")


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
        # We are older than 3.8, import from external importlib_metadata instead.
        from importlib_metadata import distributions

    pkg_to_dist = collections.defaultdict(list)
    for dist in distributions():
        for pkg in (dist.read_text('top_level.txt') or '').split():
            pkg_to_dist[pkg].append(dist.metadata['Name'])
    return dict(pkg_to_dist)
