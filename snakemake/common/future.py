"""
Future features of Python already available here.
"""
__author__ = "Hielke Walinga"
__copyright__ = "Copyright 2022, Hielke Walinga"
__email__ = "hielkewalinga@gmail.com"
__license__ = "MIT"

import collections
import importlib
from _thread import RLock


def packages_distributions():
    """Return a mapping of top-level packages to their distributions."""

    # Source: https://github.com/python/cpython/blob/9006b4471cc4d74d0cbf782d841a330b0f46a3d0/Lib/importlib/metadata/__init__.py#L1031

    # TODO: Transition instructions:
    #
    # Deprecation of 3.7
    # importlib.metadata becomes part of stdlib from 3.8 onwards
    # importlib_metadata can be removed from deps
    #
    # Deprecation of 3.9
    # importlib.metadata will have package_distributions from 3.10 onwards

    try:
        from importlib.metadata import packages_distributions as pkgs_dist

        return pkgs_dist()
    except ImportError:
        pass  # We are older than 3.10, import distributions instead.

    if hasattr(importlib, "metadata"):
        from importlib.metadata import distributions
    else:
        # We are older than 3.8,
        # import from external importlib_metadata instead.
        from importlib_metadata import distributions

    pkg_to_dist = collections.defaultdict(list)
    for dist in distributions():
        for pkg in (dist.read_text("top_level.txt") or "").split():
            pkg_to_dist[pkg].append(dist.metadata["Name"])
    return dict(pkg_to_dist)


# Source: https://github.com/python/cpython/blob/625f6704c0d783360574bbab2f78b0b9bbed5891/Lib/functools.py#L945

# TODO: Transition instructions:
#
# Deprecation of 3.7
# cached_property becomes part of functools from 3.8 onwards


################################################################################
### cached_property() - computed once per instance, cached as attribute
################################################################################

_NOT_FOUND = object()


class cached_property:
    def __init__(self, func):
        self.func = func
        self.attrname = None
        self.__doc__ = func.__doc__
        self.lock = RLock()

    def __set_name__(self, owner, name):
        if self.attrname is None:
            self.attrname = name
        elif name != self.attrname:
            raise TypeError(
                "Cannot assign the same cached_property to two different names "
                f"({self.attrname!r} and {name!r})."
            )

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        if self.attrname is None:
            raise TypeError(
                "Cannot use cached_property instance without calling __set_name__ on it.")
        try:
            cache = instance.__dict__
        except AttributeError:  # not all objects have __dict__ (e.g. class defines slots)
            msg = (
                f"No '__dict__' attribute on {type(instance).__name__!r} "
                f"instance to cache {self.attrname!r} property."
            )
            raise TypeError(msg) from None
        val = cache.get(self.attrname, _NOT_FOUND)
        if val is _NOT_FOUND:
            with self.lock:
                # check if another thread filled cache while we awaited lock
                val = cache.get(self.attrname, _NOT_FOUND)
                if val is _NOT_FOUND:
                    val = self.func(instance)
                    try:
                        cache[self.attrname] = val
                    except TypeError:
                        msg = (
                            f"The '__dict__' attribute on {type(instance).__name__!r} instance "
                            f"does not support item assignment for caching {self.attrname!r} property."
                        )
                        raise TypeError(msg) from None
        return val
