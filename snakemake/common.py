__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

from functools import update_wrapper
import inspect
import uuid
import os


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


MIN_PY_VERSION = (3, 5)
DYNAMIC_FILL = "__snakemake_dynamic__"
SNAKEMAKE_SEARCHPATH = os.path.dirname(os.path.dirname(__file__))
UUID_NAMESPACE = uuid.uuid5(uuid.NAMESPACE_URL, "https://snakemake.readthedocs.io")


def get_last_stable_version():
    return __version__.split("+")[0]


def get_container_image():
    return "quay.io/snakemake/snakemake:v{}".format(get_last_stable_version())


def get_uuid(name):
    return uuid.uuid5(UUID_NAMESPACE, name)


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


def log_location(msg):
    callerframerecord = inspect.stack()[1]
    frame = callerframerecord[0]
    info = inspect.getframeinfo(frame)
    logger.debug("{}: {info.filename}, {info.function}, {info.lineno}".format(msg, info=info))


def escape_backslash(path):
    return path.replace("\\", "\\\\")
