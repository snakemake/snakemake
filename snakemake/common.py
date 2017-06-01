__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

from functools import update_wrapper

DYNAMIC_FILL = "__snakemake_dynamic__"


class Mode:
    """
    Enum for execution mode of Snakemake.
    This handles the behavior of e.g. the logger.
    """
    default = 0
    subprocess = 1
    cluster = 2


class lazy_property(property):
    __slots__ = ["method", "cached", "__doc__"]

    def __init__(self, method):
        self.method = method
        self.cached = "_{}".format(method.__name__)
        super().__init__(method, doc=method.__doc__)

    def __get__(self, instance, owner):
        cached = getattr(instance, self.cached) if hasattr(instance, self.cached) else None
        if cached is not None:
            return cached
        value = self.method(instance)
        setattr(instance, self.cached, value)
        return value


def strip_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text
