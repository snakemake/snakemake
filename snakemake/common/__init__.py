__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

import contextlib
import itertools
import math
import operator
import platform
import hashlib
import inspect
import shutil
import sys
import uuid
import os
import asyncio
import collections
from pathlib import Path

from snakemake import __version__
from snakemake_interface_common.exceptions import WorkflowError


MIN_PY_VERSION = (3, 7)
UUID_NAMESPACE = uuid.uuid5(uuid.NAMESPACE_URL, "https://snakemake.readthedocs.io")
NOTHING_TO_BE_DONE_MSG = (
    "Nothing to be done (all requested files are present and up to date)."
)

ON_WINDOWS = platform.system() == "Windows"
# limit the number of input/output files list in job properties
# see https://github.com/snakemake/snakemake/issues/2097
IO_PROP_LIMIT = 100
SNAKEFILE_CHOICES = list(
    map(
        Path,
        (
            "Snakefile",
            "snakefile",
            "workflow/Snakefile",
            "workflow/snakefile",
        ),
    )
)


def get_snakemake_searchpaths():
    paths = [str(Path(__file__).parent.parent.parent)] + [
        path for path in sys.path if os.path.isdir(path)
    ]
    return list(unique_justseen(paths))


def mb_to_mib(mb):
    return int(math.ceil(mb * 0.95367431640625))


def parse_key_value_arg(arg, errmsg, strip_quotes=True):
    try:
        key, val = arg.split("=", 1)
    except ValueError:
        raise ValueError(errmsg + f" (Unparsable value: {repr(arg)})")
    if strip_quotes:
        val = val.strip("'\"")
    return key, val


def dict_to_key_value_args(
    some_dict: dict, quote_str: bool = True, repr_obj: bool = False
):
    items = []
    for key, value in some_dict.items():
        if repr_obj and not isinstance(value, str):
            encoded = repr(value)
        else:
            encoded = f"'{value}'" if quote_str and isinstance(value, str) else value
        items.append(f"{key}={encoded}")
    return items


def async_run(coroutine):
    """Attaches to running event loop or creates a new one to execute a
    coroutine.
    .. seealso::
         https://github.com/snakemake/snakemake/issues/1105
         https://stackoverflow.com/a/65696398
    """
    try:
        return asyncio.run(coroutine)
    except RuntimeError as e:
        coroutine.close()
        raise WorkflowError(
            "Error running coroutine in event loop. Snakemake currently does not "
            "support being executed from an already running event loop. "
            "If you run Snakemake e.g. from a Jupyter notebook, make sure to spawn a "
            "separate process for Snakemake.",
            e,
        )


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

        uri = parse_uri(f"{base}/{path}")
        if not ON_WINDOWS:
            # Norm the path such that it does not contain any ../,
            # which is invalid in an URL.
            assert uri.uri_path[0] == "/"
            uri_path = os.path.normpath(uri.uri_path)
        else:
            uri_path = uri.uri_path
        return f"{uri.scheme}:/{uri_path}"


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
    return f"snakemake/snakemake:v{get_last_stable_version()}"


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

    def __init__(self):
        self._rules = dict()

    def _register_rule(self, name, rule):
        self._rules[name] = rule

    def __getattr__(self, name):
        from snakemake.exceptions import WorkflowError

        try:
            return self._rules[name]
        except KeyError:
            raise WorkflowError(
                f"Rule {name} is not defined in this workflow. "
                f"Available rules: {', '.join(self._rules)}"
            )


class Scatter:
    """A namespace for scatter to allow items to be accessed via dot notation."""

    pass


class Gather:
    """A namespace for gather to allow items to be accessed via dot notation."""

    pass


def get_function_params(func):
    return inspect.signature(func).parameters


def get_input_function_aux_params(func, candidate_params):
    func_params = get_function_params(func)
    has_var_keyword = any(
        param.kind == param.VAR_KEYWORD for param in func_params.values()
    )
    if has_var_keyword:
        # If the function has a **kwargs parameter, we assume that it can take any
        # parameter, so we return all candidate parameters.
        return candidate_params
    else:
        return {k: v for k, v in candidate_params.items() if k in func_params}


def unique_justseen(iterable, key=None):
    """
    List unique elements, preserving order. Remember only the element just seen.

    From https://docs.python.org/3/library/itertools.html#itertools-recipes
    """
    # unique_justseen('AAAABBBCCDAABBB') --> A B C D A B
    # unique_justseen('ABBcCAD', str.lower) --> A B c A D
    return map(next, map(operator.itemgetter(1), itertools.groupby(iterable, key)))


# Taken from https://stackoverflow.com/a/34333710/7070491.
# Thanks to Laurent Laporte.
@contextlib.contextmanager
def set_env(**environ):
    """
    Temporarily set the process environment variables.

    >>> with set_env(PLUGINS_DIR='test/plugins'):
    ...   "PLUGINS_DIR" in os.environ
    True

    >>> "PLUGINS_DIR" in os.environ
    False

    :type environ: Dict[str, unicode]
    :param environ: Environment variables to set
    """
    old_environ = dict(os.environ)
    os.environ.update(environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(old_environ)


def expand_vars_and_user(value):
    if value is not None:
        return os.path.expanduser(os.path.expandvars(value))


# Taken from https://stackoverflow.com/a/2166841/7070491
# Thanks to Alex Martelli.
def is_namedtuple_instance(x):
    t = type(x)
    b = t.__bases__
    if len(b) != 1 or b[0] != tuple:
        return False
    f = getattr(t, "_fields", None)
    if not isinstance(f, tuple):
        return False
    return all(type(n) is str for n in f)


def copy_permission_safe(src: str, dst: str):
    """Copy a file to a given destination.

    If destination exists, it is removed first in order to avoid permission issues when
    the destination permissions are tried to be applied to an already existing
    destination.
    """
    if os.path.exists(dst):
        os.unlink(dst)
    shutil.copy(src, dst)
