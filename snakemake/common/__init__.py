__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

from functools import update_wrapper
import itertools
import platform
import hashlib
import inspect
import uuid
import os
import asyncio
import sys
import collections
from pathlib import Path

from snakemake._version import get_versions

__version__ = get_versions()["version"]
del get_versions


MIN_PY_VERSION = (3, 5)
DYNAMIC_FILL = "__snakemake_dynamic__"
SNAKEMAKE_SEARCHPATH = str(Path(__file__).parent.parent.parent)
UUID_NAMESPACE = uuid.uuid5(uuid.NAMESPACE_URL, "https://snakemake.readthedocs.io")
NOTHING_TO_BE_DONE_MSG = (
    "Nothing to be done (all requested files are present and up to date)."
)

ON_WINDOWS = platform.system() == "Windows"


if sys.version_info < (3, 7):

    def async_run(coroutine):
        loop = asyncio.get_event_loop()
        return loop.run_until_complete(coroutine)

else:

    def async_run(coroutine):
        """Attaches to running event loop or creates a new one to execute a
        coroutine.

        .. seealso::

             https://github.com/snakemake/snakemake/issues/1105
             https://stackoverflow.com/a/65696398

        """
        try:
            _ = asyncio.get_running_loop()
        except RuntimeError:
            asyncio.run(coroutine)
        else:
            asyncio.create_task(coroutine)


# A string that prints as TBD
class TBDString(str):
    # the second arg is necessary to avoid problems when pickling
    def __new__(cls, _=None):
        return str.__new__(cls, "<TBD>")


APPDIRS = None


RULEFUNC_CONTEXT_MARKER = "__is_snakemake_rule_func"


def get_appdirs():
    global APPDIRS
    if APPDIRS is None:
        from appdirs import AppDirs

        APPDIRS = AppDirs("snakemake", "snakemake")
    return APPDIRS


def is_local_file(path_or_uri):
    return parse_uri(path_or_uri).scheme == "file"


def parse_uri(path_or_uri):
    from smart_open import parse_uri

    try:
        return parse_uri(path_or_uri)
    except NotImplementedError as e:
        # Snakemake sees a lot of URIs which are not supported by smart_open yet
        # "docker", "git+file", "shub", "ncbi","root","roots","rootk", "gsiftp",
        # "srm","ega","ab","dropbox"
        # Fall back to a simple split if we encounter something which isn't supported.
        scheme, _, uri_path = path_or_uri.partition("://")
        if scheme and uri_path:
            uri = collections.namedtuple("Uri", ["scheme", "uri_path"])
            return uri(scheme, uri_path)
        else:
            raise e


def smart_join(base, path, abspath=False):
    if is_local_file(base):
        full = os.path.join(base, path)
        if abspath:
            return os.path.abspath(full)
        return full
    else:
        from smart_open import parse_uri

        uri = parse_uri("{}/{}".format(base, path))
        if not ON_WINDOWS:
            # Norm the path such that it does not contain any ../,
            # which is invalid in an URL.
            assert uri.uri_path[0] == "/"
            uri_path = os.path.normpath(uri.uri_path)
        else:
            uri_path = uri.uri_path
        return "{scheme}:/{uri_path}".format(scheme=uri.scheme, uri_path=uri_path)


def num_if_possible(s):
    """Convert string to number if possible, otherwise return string."""
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def get_last_stable_version():
    return __version__.split("+")[0]


def get_container_image():
    return "snakemake/snakemake:v{}".format(get_last_stable_version())


def get_uuid(name):
    return uuid.uuid5(UUID_NAMESPACE, name)


def get_file_hash(filename, algorithm="sha256"):
    """find the SHA256 hash string of a file. We use this so that the
    user can choose to cache working directories in storage.
    """
    from snakemake.logging import logger

    # The algorithm must be available
    try:
        hasher = hashlib.new(algorithm)
    except ValueError as ex:
        logger.error("%s is not an available algorithm." % algorithm)
        raise ex

    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def bytesto(bytes, to, bsize=1024):
    """convert bytes to megabytes.
    bytes to mb: bytesto(bytes, 'm')
    bytes to gb: bytesto(bytes, 'g' etc.
    From https://gist.github.com/shawnbutts/3906915
    """
    levels = {"k": 1, "m": 2, "g": 3, "t": 4, "p": 5, "e": 6}
    answer = float(bytes)
    for _ in range(levels[to]):
        answer = answer / bsize
    return answer


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

    @staticmethod
    def clean(instance, method):
        delattr(instance, method)

    def __init__(self, method):
        self.method = method
        self.cached = "_{}".format(method.__name__)
        super().__init__(method, doc=method.__doc__)

    def __get__(self, instance, owner):
        cached = (
            getattr(instance, self.cached) if hasattr(instance, self.cached) else None
        )
        if cached is not None:
            return cached
        value = self.method(instance)
        setattr(instance, self.cached, value)
        return value


def strip_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix) :]
    return text


def log_location(msg):
    from snakemake.logging import logger

    callerframerecord = inspect.stack()[1]
    frame = callerframerecord[0]
    info = inspect.getframeinfo(frame)
    logger.debug(
        "{}: {info.filename}, {info.function}, {info.lineno}".format(msg, info=info)
    )


def group_into_chunks(n, iterable):
    """Group iterable into chunks of size at most n.

    See https://stackoverflow.com/a/8998040.
    """
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk


class Rules:
    """A namespace for rules so that they can be accessed via dot notation."""

    pass


class Scatter:
    """A namespace for scatter to allow items to be accessed via dot notation."""

    pass


class Gather:
    """A namespace for gather to allow items to be accessed via dot notation."""

    pass
