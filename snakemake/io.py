__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import collections
import os
import shutil
from pathlib import Path
import re
import stat
import time
import datetime
import json
import copy
import functools
import subprocess as sp
from itertools import product, chain
from contextlib import contextmanager
import string
import collections
import asyncio

from snakemake.exceptions import (
    MissingOutputException,
    WorkflowError,
    WildcardError,
    RemoteFileException,
)
from snakemake.logging import logger
from inspect import isfunction, ismethod
from snakemake.common import DYNAMIC_FILL, ON_WINDOWS, async_run


class Mtime:
    __slots__ = ["_local", "_local_target", "_remote"]

    def __init__(self, local=None, local_target=None, remote=None):
        self._local = local
        self._local_target = local_target
        self._remote = remote

    def local_or_remote(self, follow_symlinks=False):
        if self._remote is not None:
            return self._remote
        if follow_symlinks and self._local_target is not None:
            return self._local_target
        return self._local

    def remote(
        self,
    ):
        return self._remote

    def local(self, follow_symlinks=False):
        if follow_symlinks and self._local_target is not None:
            return self._local_target
        return self._local


def lutime(f, times):
    # In some cases, we have a platform where os.supports_follow_symlink includes stat()
    # but not utime().  This leads to an anomaly.  In any case we never want to touch the
    # target of a link.
    if os.utime in os.supports_follow_symlinks:
        # ...utime is well behaved
        os.utime(f, times, follow_symlinks=False)
    elif not os.path.islink(f):
        # ...symlinks not an issue here
        os.utime(f, times)
    else:
        try:
            # try the system command
            if times:
                fmt_time = lambda sec: datetime.fromtimestamp(sec).strftime(
                    "%Y%m%d%H%M.%S"
                )
                atime, mtime = times
                sp.check_call(["touch", "-h", f, "-a", "-t", fmt_time(atime)])
                sp.check_call(["touch", "-h", f, "-m", "-t", fmt_time(mtime)])
            else:
                sp.check_call(["touch", "-h", f])
        except sp.CalledProcessError:
            pass
        # ...problem system.  Do nothing.
        logger.warning(
            "Unable to set utime on symlink {}. Your Python build does not support it.".format(
                f
            )
        )
        return None


if os.chmod in os.supports_follow_symlinks:

    def lchmod(f, mode):
        os.chmod(f, mode, follow_symlinks=False)


else:

    def lchmod(f, mode):
        os.chmod(f, mode)


class ExistsDict(dict):
    def __init__(self, cache):
        super().__init__()
        self.cache = cache
        self.has_inventory = set()

    def __getitem__(self, path):
        # Always return False if not in dict.
        # The reason is that this is only called if the method contains below has returned True.
        # Hence, we already know that either path is in dict, or inventory has never
        # seen it, and hence it does not exist.
        return self.get(path, False)

    def __contains__(self, path):
        # if already in inventory, always return True.
        parent = path.get_inventory_parent()
        return parent in self.has_inventory or super().__contains__(path)


class IOCache:
    def __init__(self, max_wait_time):
        self.mtime = dict()
        self.exists_local = ExistsDict(self)
        self.exists_remote = ExistsDict(self)
        self.size = dict()
        self.active = True
        self.remaining_wait_time = max_wait_time
        self.max_wait_time = max_wait_time

    def mtime_inventory(self, jobs):
        async_run(self._mtime_inventory(jobs))

    async def _mtime_inventory(self, jobs, n_workers=8):
        queue = asyncio.Queue()
        stop_item = object()

        async def worker(queue):
            while True:
                item = await queue.get()
                if item is stop_item:
                    queue.task_done()
                    return
                try:
                    self.mtime[item] = await self.collect_mtime(item)
                except Exception as e:
                    queue.task_done()

                    raise e
                queue.task_done()

        tasks = [
            asyncio.get_event_loop().create_task(worker(queue))
            for _ in range(n_workers)
        ]

        for job in jobs:
            for f in chain(job.input, job.expanded_output):
                if f.exists:
                    queue.put_nowait(f)
            if job.benchmark and job.benchmark.exists:
                queue.put_nowait(job.benchmark)

        # Send a stop item to each worker.
        for _ in range(n_workers):
            queue.put_nowait(stop_item)

        await asyncio.gather(*tasks)

    async def collect_mtime(self, path):
        return path.mtime_uncached

    def clear(self):
        self.mtime.clear()
        self.size.clear()
        self.exists_local.clear()
        self.exists_remote.clear()
        self.remaining_wait_time = self.max_wait_time

    def deactivate(self):
        self.clear()
        self.active = False


def IOFile(file, rule=None):
    assert rule is not None
    f = _IOFile(file)
    f.rule = rule
    return f


class _IOFile(str):
    """
    A file that is either input or output of a rule.
    """

    __slots__ = [
        "_is_function",
        "_file",
        "rule",
        "_regex",
    ]

    def __new__(cls, file):
        is_annotated = isinstance(file, AnnotatedString)
        is_callable = (
            isfunction(file) or ismethod(file) or (is_annotated and bool(file.callable))
        )
        if not is_callable and file.endswith("/"):
            # remove trailing slashes
            stripped = file.rstrip("/")
            if is_annotated:
                stripped = AnnotatedString(stripped)
                stripped.flags = file.flags
            file = stripped
        obj = str.__new__(cls, file)
        obj._is_function = is_callable
        obj._file = file
        obj.rule = None
        obj._regex = None

        if obj.is_remote:
            obj.remote_object._iofile = obj

        return obj

    def iocache(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            if self.rule.workflow.iocache.active:
                cache = getattr(self.rule.workflow.iocache, func.__name__)
                if self in cache:
                    return cache[self]
                v = func(self, *args, **kwargs)
                cache[self] = v
                return v
            else:
                return func(self, *args, **kwargs)

        return wrapper

    def _refer_to_remote(func):
        """
        A decorator so that if the file is remote and has a version
        of the same file-related function, call that version instead.
        """

        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            if self.is_remote:
                if hasattr(self.remote_object, func.__name__):
                    return getattr(self.remote_object, func.__name__)(*args, **kwargs)
            return func(self, *args, **kwargs)

        return wrapper

    def inventory(self):
        async_run(self._inventory())

    async def _inventory(self):
        """Starting from the given file, try to cache as much existence and
        modification date information of this and other files as possible.
        """
        cache = self.rule.workflow.iocache
        if cache.active:
            tasks = []
            if self.is_remote and self not in cache.exists_remote:
                # info not yet in inventory, let's discover as much as we can
                tasks.append(self.remote_object.inventory(cache))
            if not ON_WINDOWS and self not in cache.exists_local:
                # we don't want to mess with different path representations on windows
                tasks.append(self._local_inventory(cache))
            await asyncio.gather(*tasks)

    async def _local_inventory(self, cache):
        # for local files, perform BFS via os.scandir to determine existence of files
        if cache.remaining_wait_time <= 0:
            # No more time to create inventory.
            return

        start_time = time.time()

        folders = self.split("/")[:-1]
        if not folders:
            return

        if os.path.isabs(self):
            # For absolute paths, only use scan the immediate parent
            ancestors = [os.path.dirname(self)]
        else:
            ancestors = ["/".join(folders[:i]) for i in range(1, len(folders) + 1)]

        for (i, path) in enumerate(ancestors):
            if path in cache.exists_local.has_inventory:
                # This path was already scanned before, hence we can stop.
                break
            try:
                with os.scandir(path) as scan:
                    for entry in scan:
                        cache.exists_local[entry.path] = True
                cache.exists_local[path] = True
                cache.exists_local.has_inventory.add(path)
            except FileNotFoundError:
                # Not found, hence, all subfolders cannot be present as well
                for path in ancestors[i:]:
                    cache.exists_local[path] = False
                    cache.exists_local.has_inventory.add(path)
                break
            except PermissionError:
                raise WorkflowError(
                    "Insufficient permissions to access {}. "
                    "Please make sure that all accessed files and directories "
                    "are readable and writable for you.".format(self)
                )

        cache.remaining_wait_time -= time.time() - start_time

    @_refer_to_remote
    def get_inventory_parent(self):
        """If eligible for inventory, get the parent of a given path.

        This code does not work on local Windows paths,
        but inventory is disabled on Windows.
        """
        parent = os.path.dirname(self)
        if parent and parent != "..":
            return parent

    @contextmanager
    def open(self, mode="r", buffering=-1, encoding=None, errors=None, newline=None):
        """Open this file. If necessary, download it from remote first.

        This can (and should) be used in a `with`-statement.
        """
        if not self.exists:
            raise WorkflowError(
                "File {} cannot be opened, since it does not exist.".format(self)
            )
        if not self.exists_local and self.is_remote:
            self.download_from_remote()

        f = open(self)
        try:
            yield f
        finally:
            f.close()

    def contains_wildcard(self):
        return contains_wildcard(self.file)

    @property
    def is_remote(self):
        return is_flagged(self._file, "remote_object")

    @property
    def is_ancient(self):
        return is_flagged(self._file, "ancient")

    @property
    def is_directory(self):
        return is_flagged(self._file, "directory")

    @property
    def is_temp(self):
        return is_flagged(self._file, "temp")

    @property
    def is_multiext(self):
        return is_flagged(self._file, "multiext")

    @property
    def multiext_prefix(self):
        return get_flag_value(self._file, "multiext")

    def update_remote_filepath(self):
        # if the file string is different in the iofile, update the remote object
        # (as in the case of wildcard expansion)
        remote_object = self.remote_object
        if remote_object._file != self._file:
            remote_object._iofile = self

    @property
    def should_keep_local(self):
        return self.remote_object.keep_local

    @property
    def should_stay_on_remote(self):
        return self.remote_object.stay_on_remote

    @property
    def remote_object(self):
        return get_flag_value(self._file, "remote_object")

    @property
    @_refer_to_remote
    def file(self):
        if not self._is_function:
            return self._file
        else:
            raise ValueError(
                "This IOFile is specified as a function and "
                "may not be used directly."
            )

    def check(self):
        hint = (
            "It can also lead to inconsistent results of the file-matching "
            "approach used by Snakemake."
        )
        if self._file.startswith("./"):
            logger.warning(
                "Relative file path '{}' starts with './'. This is redundant "
                "and strongly discouraged. {} You can simply omit the './' "
                "for relative file paths.".format(self._file, hint)
            )
        if self._file.startswith(" "):
            logger.warning(
                "File path '{}' starts with whitespace. "
                "This is likely unintended. {}".format(self._file, hint)
            )
        if self._file.endswith(" "):
            logger.warning(
                "File path '{}' ends with whitespace. "
                "This is likely unintended. {}".format(self._file, hint)
            )
        if "\n" in self._file:
            logger.warning(
                "File path '{}' contains line break. "
                "This is likely unintended. {}".format(self._file, hint)
            )
        if _double_slash_regex.search(self._file) is not None and not self.is_remote:
            logger.warning(
                "File path {} contains double '{}'. "
                "This is likely unintended. {}".format(self._file, os.path.sep, hint)
            )

    @property
    def exists(self):
        if self.is_remote:
            return self.exists_remote
        else:
            return self.exists_local

    def parents(self, omit=0):
        """Yield all parent paths, omitting the given number of ancestors."""
        for p in list(Path(self.file).parents)[::-1][omit:]:
            p = IOFile(str(p), rule=self.rule)
            p.clone_flags(self)
            yield p

    @property
    @iocache
    def exists_local(self):
        return os.path.exists(self.file)

    @property
    @iocache
    def exists_remote(self):
        if not self.is_remote:
            return False
        return self.remote_object.exists()

    @property
    def protected(self):
        """Returns True if the file is protected. Always False for symlinks."""
        # symlinks are never regarded as protected
        return (
            self.exists_local
            and not os.access(self.file, os.W_OK)
            and not os.path.islink(self.file)
        )

    @property
    @iocache
    def mtime(self):
        return self.mtime_uncached

    @property
    def mtime_uncached(self):
        """Obtain mtime.

        Usually, this will be one stat call only. For symlinks and directories
        it will be two, for symlinked directories it will be three,
        for remote files it will additionally query the remote
        location.
        """
        mtime_remote = self.remote_object.mtime() if self.is_remote else None

        # We first do a normal stat.
        try:
            _stat = os.stat(self.file, follow_symlinks=False)

            is_symlink = stat.S_ISLNK(_stat.st_mode)
            is_dir = stat.S_ISDIR(_stat.st_mode)
            mtime = _stat.st_mtime

            def get_dir_mtime():
                # Try whether we have a timestamp file for it.
                return os.stat(
                    os.path.join(self.file, ".snakemake_timestamp"),
                    follow_symlinks=True,
                ).st_mtime

            if not is_symlink:
                if is_dir:
                    try:
                        mtime = get_dir_mtime()
                    except FileNotFoundError:
                        # No timestamp, hence go on as if it is a file.
                        pass

                # In the usual case, not a dir, not a symlink.
                # We need just a single stat call.
                return Mtime(local=mtime, remote=mtime_remote)

            else:
                # In case of a symlink, we need the stats for the target file/dir.
                target_stat = os.stat(self.file, follow_symlinks=True)
                # Further, we need to check again if this is a directory.
                is_dir = stat.S_ISDIR(target_stat.st_mode)
                mtime_target = target_stat.st_mtime

                if is_dir:
                    try:
                        mtime_target = get_dir_mtime()
                    except FileNotFoundError:
                        # No timestamp, hence go on as if it is a file.
                        pass

                return Mtime(
                    local=mtime, local_target=mtime_target, remote=mtime_remote
                )

        except FileNotFoundError:
            if self.is_remote:
                return Mtime(remote=mtime_remote)
            raise WorkflowError(
                "Unable to obtain modification time of file {} although it existed before. "
                "It could be that a concurrent process has deleted it while Snakemake "
                "was running.".format(self.file)
            )
        except PermissionError:
            raise WorkflowError(
                "Unable to obtain modification time of file {} because of missing "
                "read permissions.".format(self.file)
            )

    @property
    def flags(self):
        return getattr(self._file, "flags", {})

    @property
    @iocache
    @_refer_to_remote
    def size(self):
        return self.size_local

    @property
    def size_local(self):
        # follow symlinks but throw error if invalid
        self.check_broken_symlink()
        return os.path.getsize(self.file)

    def check_broken_symlink(self):
        """ Raise WorkflowError if file is a broken symlink. """
        if not self.exists_local and os.lstat(self.file):
            raise WorkflowError(
                "File {} seems to be a broken symlink.".format(self.file)
            )

    @_refer_to_remote
    def is_newer(self, time):
        """Returns true of the file (which is an input file) is newer than time, or if it is
        a symlink that points to a file newer than time."""
        if self.is_ancient:
            return False

        return self.mtime.local_or_remote(follow_symlinks=True) > time

    def download_from_remote(self):
        if self.is_remote and self.remote_object.exists():
            if not self.should_stay_on_remote:
                logger.info("Downloading from remote: {}".format(self.file))
                self.remote_object.download()
                logger.info("Finished download.")
        else:
            raise RemoteFileException(
                "The file to be downloaded does not seem to exist remotely."
            )

    def upload_to_remote(self):
        if self.is_remote:
            logger.info("Uploading to remote: {}".format(self.file))
            self.remote_object.upload()
            logger.info("Finished upload.")

    def prepare(self):
        path_until_wildcard = re.split(DYNAMIC_FILL, self.file)[0]
        dir = os.path.dirname(path_until_wildcard)
        if len(dir) > 0:
            try:
                os.makedirs(dir, exist_ok=True)
            except OSError as e:
                # ignore Errno 17 "File exists" (reason: multiprocessing)
                if e.errno != 17:
                    raise e

        if is_flagged(self._file, "pipe"):
            os.mkfifo(self._file)

    def protect(self):
        mode = (
            os.lstat(self.file).st_mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH
        )
        if os.path.isdir(self.file):
            for root, dirs, files in os.walk(self.file):
                for d in dirs:
                    lchmod(os.path.join(self.file, d), mode)
                for f in files:
                    lchmod(os.path.join(self.file, f), mode)
        lchmod(self.file, mode)

    def remove(self, remove_non_empty_dir=False):
        if self.is_directory:
            remove(self, remove_non_empty_dir=True)
        else:
            remove(self, remove_non_empty_dir=remove_non_empty_dir)

    def touch(self, times=None):
        """ times must be 2-tuple: (atime, mtime) """
        try:
            if self.is_directory:
                file = os.path.join(self.file, ".snakemake_timestamp")
                # Create the flag file if it doesn't exist
                if not os.path.exists(file):
                    with open(file, "w"):
                        pass
                lutime(file, times)
            else:
                lutime(self.file, times)
        except OSError as e:
            if e.errno == 2:
                raise MissingOutputException(
                    "Output file {} of rule {} shall be touched but "
                    "does not exist.".format(self.file, self.rule.name),
                    lineno=self.rule.lineno,
                    snakefile=self.rule.snakefile,
                )
            else:
                raise e

    def touch_or_create(self):
        try:
            self.touch()
        except MissingOutputException:
            # first create directory if it does not yet exist
            dir = self.file if self.is_directory else os.path.dirname(self.file)
            if dir:
                os.makedirs(dir, exist_ok=True)
            # create empty file
            file = (
                os.path.join(self.file, ".snakemake_timestamp")
                if self.is_directory
                else self.file
            )
            with open(file, "w") as f:
                pass

    def apply_wildcards(self, wildcards, fill_missing=False, fail_dynamic=False):
        f = self._file

        if self._is_function:
            f = self._file(Namedlist(fromdict=wildcards))

        # this bit ensures flags are transferred over to files after
        # wildcards are applied

        file_with_wildcards_applied = IOFile(
            apply_wildcards(
                f,
                wildcards,
                fill_missing=fill_missing,
                fail_dynamic=fail_dynamic,
                dynamic_fill=DYNAMIC_FILL,
            ),
            rule=self.rule,
        )

        file_with_wildcards_applied.clone_flags(self)

        return file_with_wildcards_applied

    def get_wildcard_names(self):
        return get_wildcard_names(self.file)

    def regex(self):
        if self._regex is None:
            # compile a regular expression
            self._regex = re.compile(regex(self.file))
        return self._regex

    def constant_prefix(self):
        first_wildcard = _wildcard_regex.search(self.file)
        if first_wildcard:
            return self.file[: first_wildcard.start()]
        return self.file

    def constant_suffix(self):
        m = None
        for m in _wildcard_regex.finditer(self.file):
            pass
        last_wildcard = m
        if last_wildcard:
            return self.file[last_wildcard.end() :]
        return self.file

    def match(self, target):
        return self.regex().match(target) or None

    def format_dynamic(self):
        return self.replace(DYNAMIC_FILL, "{*}")

    def clone_flags(self, other):
        if isinstance(self._file, str):
            self._file = AnnotatedString(self._file)
        if isinstance(other._file, AnnotatedString):
            self._file.flags = getattr(other._file, "flags", {}).copy()
            if "remote_object" in self._file.flags:
                self._file.flags["remote_object"] = copy.copy(
                    self._file.flags["remote_object"]
                )
                self.update_remote_filepath()

    def clone_remote_object(self, other):
        if (
            isinstance(other._file, AnnotatedString)
            and "remote_object" in other._file.flags
        ):
            self._file.flags["remote_object"] = copy.copy(
                other._file.flags["remote_object"]
            )
            self.update_remote_filepath()

    def set_flags(self, flags):
        if isinstance(self._file, str):
            self._file = AnnotatedString(self._file)
        self._file.flags = flags

    def __eq__(self, other):
        f = other._file if isinstance(other, _IOFile) else other
        return self._file == f

    def __hash__(self):
        return self._file.__hash__()


_double_slash_regex = (
    re.compile(r"([^:]//|^//)") if os.path.sep == "/" else re.compile(r"\\\\")
)


_wildcard_regex = re.compile(
    r"""
    \{
        (?=(   # This lookahead assertion emulates an 'atomic group'
               # which is required for performance
            \s*(?P<name>\w+)                    # wildcard name
            (\s*,\s*
                (?P<constraint>                 # an optional constraint
                    ([^{}]+ | \{\d+(,\d+)?\})*  # allow curly braces to nest one level
                )                               # ...  as in '{w,a{3,5}}'
            )?\s*
        ))\1
    \}
    """,
    re.VERBOSE,
)


def wait_for_files(
    files, latency_wait=3, force_stay_on_remote=False, ignore_pipe=False
):
    """Wait for given files to be present in filesystem."""
    files = list(files)

    def get_missing():
        return [
            f
            for f in files
            if not (
                f.exists_remote
                if (
                    isinstance(f, _IOFile)
                    and f.is_remote
                    and (force_stay_on_remote or f.should_stay_on_remote)
                )
                else os.path.exists(f)
                if not (is_flagged(f, "pipe") and ignore_pipe)
                else True
            )
        ]

    missing = get_missing()
    if missing:
        logger.info(
            "Waiting at most {} seconds for missing files.".format(latency_wait)
        )
        for _ in range(latency_wait):
            if not get_missing():
                return
            time.sleep(1)
        raise IOError(
            "Missing files after {} seconds:\n{}".format(
                latency_wait, "\n".join(get_missing())
            )
        )


def get_wildcard_names(pattern):
    return set(match.group("name") for match in _wildcard_regex.finditer(pattern))


def contains_wildcard(path):
    return _wildcard_regex.search(path) is not None


def contains_wildcard_constraints(pattern):
    return any(match.group("constraint") for match in _wildcard_regex.finditer(pattern))


def remove(file, remove_non_empty_dir=False):
    if file.is_remote and file.should_stay_on_remote:
        if file.exists_remote:
            file.remote_object.remove()
    elif os.path.isdir(file) and not os.path.islink(file):
        if remove_non_empty_dir:
            shutil.rmtree(file)
        else:
            try:
                os.removedirs(file)
            except OSError as e:
                # skip non empty directories
                if e.errno == 39:
                    logger.info(
                        "Skipped removing non-empty directory {}".format(e.filename)
                    )
                else:
                    logger.warning(str(e))
    # Remember that dangling symlinks fail the os.path.exists() test, but
    # we definitely still want to zap them. try/except is the safest way.
    # Also, we don't want to remove the null device if it is an output.
    elif os.devnull != str(file):
        try:
            os.remove(file)
        except FileNotFoundError:
            pass


def regex(filepattern):
    f = []
    last = 0
    wildcards = set()
    for match in _wildcard_regex.finditer(filepattern):
        f.append(re.escape(filepattern[last : match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "Constraint regex must be defined only in the first "
                    "occurence of the wildcard in a string."
                )
            f.append("(?P={})".format(wildcard))
        else:
            wildcards.add(wildcard)
            f.append(
                "(?P<{}>{})".format(
                    wildcard,
                    match.group("constraint") if match.group("constraint") else ".+",
                )
            )
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$")  # ensure that the match spans the whole file
    return "".join(f)


def apply_wildcards(
    pattern,
    wildcards,
    fill_missing=False,
    fail_dynamic=False,
    dynamic_fill=None,
    keep_dynamic=False,
):
    def format_match(match):
        name = match.group("name")
        try:
            value = wildcards[name]
            if fail_dynamic and value == dynamic_fill:
                raise WildcardError(name)
            return str(value)  # convert anything into a str
        except KeyError as ex:
            if keep_dynamic:
                return "{{{}}}".format(name)
            elif fill_missing:
                return dynamic_fill
            else:
                raise WildcardError(str(ex))

    return _wildcard_regex.sub(format_match, pattern)


def not_iterable(value):
    return (
        isinstance(value, str)
        or isinstance(value, dict)
        or not isinstance(value, collections.abc.Iterable)
    )


def is_callable(value):
    return (
        callable(value)
        or (isinstance(value, _IOFile) and value._is_function)
        or (isinstance(value, AnnotatedString) and value.callable is not None)
    )


class AnnotatedString(str):
    def __init__(self, value):
        self.flags = dict()
        self.callable = value if is_callable(value) else None


def flag(value, flag_type, flag_value=True):
    if isinstance(value, AnnotatedString):
        value.flags[flag_type] = flag_value
        return value
    if not_iterable(value):
        value = AnnotatedString(value)
        value.flags[flag_type] = flag_value
        return value
    return [flag(v, flag_type, flag_value=flag_value) for v in value]


def is_flagged(value, flag):
    if isinstance(value, AnnotatedString):
        return flag in value.flags and value.flags[flag]
    if isinstance(value, _IOFile):
        return flag in value.flags and value.flags[flag]
    return False


def get_flag_value(value, flag_type):
    if isinstance(value, AnnotatedString) or isinstance(value, _IOFile):
        if flag_type in value.flags:
            return value.flags[flag_type]
        else:
            return None


def ancient(value):
    """
    A flag for an input file that shall be considered ancient; i.e. its timestamp shall have no effect on which jobs to run.
    """
    return flag(value, "ancient")


def directory(value):
    """
    A flag to specify that an output is a directory, rather than a file or named pipe.
    """
    if is_flagged(value, "pipe"):
        raise SyntaxError("Pipe and directory flags are mutually exclusive.")
    if is_flagged(value, "remote"):
        raise SyntaxError("Remote and directory flags are mutually exclusive.")
    if is_flagged(value, "dynamic"):
        raise SyntaxError("Dynamic and directory flags are mutually exclusive.")
    return flag(value, "directory")


def temp(value):
    """
    A flag for an input or output file that shall be removed after usage.
    """
    if is_flagged(value, "protected"):
        raise SyntaxError("Protected and temporary flags are mutually exclusive.")
    if is_flagged(value, "remote"):
        raise SyntaxError("Remote and temporary flags are mutually exclusive.")
    return flag(value, "temp")


def pipe(value):
    if is_flagged(value, "protected"):
        raise SyntaxError("Pipes may not be protected.")
    if is_flagged(value, "remote"):
        raise SyntaxError("Pipes may not be remote files.")
    if ON_WINDOWS:
        logger.warning("Pipes is not yet supported on Windows.")
    return flag(value, "pipe", not ON_WINDOWS)


def temporary(value):
    """ An alias for temp. """
    return temp(value)


def protected(value):
    """ A flag for a file that shall be write protected after creation. """
    if is_flagged(value, "temp"):
        raise SyntaxError("Protected and temporary flags are mutually exclusive.")
    if is_flagged(value, "remote"):
        raise SyntaxError("Remote and protected flags are mutually exclusive.")
    return flag(value, "protected")


def dynamic(value):
    """
    A flag for a file that shall be dynamic, i.e. the multiplicity
    (and wildcard values) will be expanded after a certain
    rule has been run"""
    annotated = flag(value, "dynamic", True)
    tocheck = [annotated] if not_iterable(annotated) else annotated
    for file in tocheck:
        matches = list(_wildcard_regex.finditer(file))
        # if len(matches) != 1:
        #    raise SyntaxError("Dynamic files need exactly one wildcard.")
        for match in matches:
            if match.group("constraint"):
                raise SyntaxError(
                    "The wildcards in dynamic files cannot be constrained."
                )
    return annotated


def touch(value):
    return flag(value, "touch")


def unpack(value):
    return flag(value, "unpack")


def repeat(value, n_repeat):
    """Flag benchmark records with the number of repeats."""
    return flag(value, "repeat", n_repeat)


def checkpoint_target(value):
    return flag(value, "checkpoint_target")


ReportObject = collections.namedtuple(
    "ReportObject", ["caption", "category", "subcategory", "patterns", "htmlindex"]
)


def report(
    value, caption=None, category=None, subcategory=None, patterns=[], htmlindex=None
):
    """Flag output file or directory as to be included into reports.

    In case of directory, files to include can be specified via a glob pattern (default: *).

    Arguments
    value -- File or directory.
    caption -- Path to a .rst file with a textual description of the result.
    category -- Name of the category in which the result should be displayed in the report.
    pattern -- Wildcard pattern for selecting files if a directory is given (this is used as
               input for snakemake.io.glob_wildcards). Pattern shall not include the path to the
               directory itself.
    """
    return flag(
        value,
        "report",
        ReportObject(caption, category, subcategory, patterns, htmlindex),
    )


def local(value):
    """Mark a file as local file. This disables application of a default remote
    provider.
    """
    if is_flagged(value, "remote"):
        raise SyntaxError("Remote and local flags are mutually exclusive.")
    return flag(value, "local")


def expand(*args, **wildcards):
    """
    Expand wildcards in given filepatterns.

    Arguments
    *args -- first arg: filepatterns as list or one single filepattern,
        second arg (optional): a function to combine wildcard values
        (itertools.product per default)
    **wildcards -- the wildcards as keyword arguments
        with their values as lists. If allow_missing=True is included
        wildcards in filepattern without values will stay unformatted.
    """
    filepatterns = args[0]
    if len(args) == 1:
        combinator = product
    elif len(args) == 2:
        combinator = args[1]
    if isinstance(filepatterns, str) or isinstance(filepatterns, Path):
        filepatterns = [filepatterns]

    def path_to_str(f):
        if isinstance(f, Path):
            return str(f)
        return f

    filepatterns = list(map(path_to_str, filepatterns))

    if any(map(lambda f: getattr(f, "flags", {}), filepatterns)):
        raise WorkflowError(
            "Flags in file patterns given to expand() are invalid. "
            "Flags (e.g. temp(), directory()) have to be applied outside "
            "of expand (e.g. 'temp(expand(\"plots/{sample}.pdf\", sample=SAMPLES))')."
        )

    # check if remove missing is provided
    format_dict = dict
    if "allow_missing" in wildcards and wildcards["allow_missing"] is True:

        class FormatDict(dict):
            def __missing__(self, key):
                return "{" + key + "}"

        format_dict = FormatDict
        # check that remove missing is not a wildcard in the filepatterns
        for filepattern in filepatterns:
            if "allow_missing" in re.findall(r"{([^}\.[!:]+)", filepattern):
                format_dict = dict
                break

    # remove unused wildcards to avoid duplicate filepatterns
    wildcards = {
        filepattern: {
            k: v
            for k, v in wildcards.items()
            if k in re.findall(r"{([^}\.[!:]+)", filepattern)
        }
        for filepattern in filepatterns
    }

    def flatten(wildcards):
        for wildcard, values in wildcards.items():
            if isinstance(values, str) or not isinstance(
                values, collections.abc.Iterable
            ):
                values = [values]
            yield [(wildcard, value) for value in values]

    formatter = string.Formatter()
    try:
        return [
            formatter.vformat(filepattern, (), comb)
            for filepattern in filepatterns
            for comb in map(format_dict, combinator(*flatten(wildcards[filepattern])))
        ]
    except KeyError as e:
        raise WildcardError("No values given for wildcard {}.".format(e))


def multiext(prefix, *extensions):
    """Expand a given prefix with multiple extensions (e.g. .txt, .csv, _peaks.bed, ...)."""
    if any((r"/" in ext or r"\\" in ext) for ext in extensions):
        raise WorkflowError(
            r"Extensions for multiext may not contain path delimiters " r"(/,\)."
        )
    return [flag(prefix + ext, "multiext", flag_value=prefix) for ext in extensions]


def limit(pattern, **wildcards):
    """
    Limit wildcards to the given values.

    Arguments:
    **wildcards -- the wildcards as keyword arguments
                   with their values as lists
    """
    return pattern.format(
        **{
            wildcard: "{{{},{}}}".format(wildcard, "|".join(values))
            for wildcard, values in wildcards.items()
        }
    )


def glob_wildcards(pattern, files=None, followlinks=False):
    """
    Glob the values of the wildcards by matching the given pattern to the filesystem.
    Returns a named tuple with a list of values for each wildcard.
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    dirname = (
        os.path.dirname(pattern[: first_wildcard.start()])
        if first_wildcard
        else os.path.dirname(pattern)
    )
    if not dirname:
        dirname = "."

    names = [match.group("name") for match in _wildcard_regex.finditer(pattern)]
    Wildcards = collections.namedtuple("Wildcards", names)
    wildcards = Wildcards(*[list() for name in names])

    pattern = re.compile(regex(pattern))

    if files is None:
        files = (
            os.path.normpath(os.path.join(dirpath, f))
            for dirpath, dirnames, filenames in os.walk(
                dirname, followlinks=followlinks
            )
            for f in chain(filenames, dirnames)
        )

    for f in files:
        match = re.match(pattern, f)
        if match:
            for name, value in match.groupdict().items():
                getattr(wildcards, name).append(value)
    return wildcards


def update_wildcard_constraints(
    pattern, wildcard_constraints, global_wildcard_constraints
):
    """Update wildcard constraints

    Args:
      pattern (str): pattern on which to update constraints
      wildcard_constraints (dict): dictionary of wildcard:constraint key-value pairs
      global_wildcard_constraints (dict): dictionary of wildcard:constraint key-value pairs
    """

    def replace_constraint(match):
        name = match.group("name")
        constraint = match.group("constraint")
        newconstraint = wildcard_constraints.get(
            name, global_wildcard_constraints.get(name)
        )
        if name in examined_names:
            return match.group(0)
        examined_names.add(name)
        # Don't override if constraint already set
        if constraint is not None:
            return match.group(0)
        # Only update if a new constraint has actually been set
        elif newconstraint is not None:
            return "{{{},{}}}".format(name, newconstraint)
        else:
            return match.group(0)

    examined_names = set()
    updated = _wildcard_regex.sub(replace_constraint, pattern)

    # inherit flags
    if isinstance(pattern, AnnotatedString):
        updated = AnnotatedString(updated)
        updated.flags = dict(pattern.flags)
    return updated


def split_git_path(path):
    file_sub = re.sub(r"^git\+file:/+", "/", path)
    (file_path, version) = file_sub.split("@")
    file_path = os.path.realpath(file_path)
    root_path = get_git_root(file_path)
    if file_path.startswith(root_path):
        file_path = file_path[len(root_path) :].lstrip("/")
    return (root_path, file_path, version)


def get_git_root(path):
    """
    Args:
        path: (str) Path a to a directory/file that is located inside the repo
    Returns:
        path to root folder for git repo
    """
    import git

    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        return git_repo.git.rev_parse("--show-toplevel")
    except git.exc.NoSuchPathError:
        tail, head = os.path.split(path)
        return get_git_root_parent_directory(tail, path)


def get_git_root_parent_directory(path, input_path):
    """
    This function will recursively go through parent directories until a git
    repository is found or until no parent directories are left, in which case
    a error will be raised. This is needed when providing a path to a
    file/folder that is located on a branch/tag no currently checked out.

    Args:
        path: (str) Path a to a directory that is located inside the repo
        input_path: (str) origin path, used when raising WorkflowError
    Returns:
        path to root folder for git repo
    """
    import git

    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        return git_repo.git.rev_parse("--show-toplevel")
    except git.exc.NoSuchPathError:
        tail, head = os.path.split(path)
        if tail is None:
            raise WorkflowError(
                "Neither provided git path ({}) ".format(input_path)
                + "or parent directories contain a valid git repo."
            )
        else:
            return get_git_root_parent_directory(tail, input_path)


def git_content(git_file):
    """
    This function will extract a file from a git repository, one located on
    the filesystem.
    Expected format is git+file:///path/to/your/repo/path_to_file@@version

    Args:
      env_file (str): consist of path to repo, @, version and file information
                      Ex: git+file:////home/smeds/snakemake-wrappers/bio/fastqc/wrapper.py@0.19.3
    Returns:
        file content or None if the expected format isn't meet
    """
    import git

    if git_file.startswith("git+file:"):
        (root_path, file_path, version) = split_git_path(git_file)
        return git.Repo(root_path).git.show("{}:{}".format(version, file_path))
    else:
        raise WorkflowError(
            "Provided git path ({}) doesn't meet the "
            "expected format:".format(git_file) + ", expected format is "
            "git+file://PATH_TO_REPO/PATH_TO_FILE_INSIDE_REPO@VERSION"
        )


def strip_wildcard_constraints(pattern):
    """Return a string that does not contain any wildcard constraints."""

    def strip_constraint(match):
        return "{{{}}}".format(match.group("name"))

    return _wildcard_regex.sub(strip_constraint, pattern)


class Namedlist(list):
    """
    A list that additionally provides functions to name items. Further,
    it is hashable, however the hash does not consider the item names.
    """

    def __init__(
        self,
        toclone=None,
        fromdict=None,
        plainstr=False,
        strip_constraints=False,
        custom_map=None,
    ):
        """
        Create the object.

        Arguments
        toclone  -- another Namedlist that shall be cloned
        fromdict -- a dict that shall be converted to a
            Namedlist (keys become names)
        """
        list.__init__(self)
        self._names = dict()

        # white-list of attribute names that can be overridden in _set_name
        # default to throwing exception if called to prevent use as functions
        self._allowed_overrides = ["index", "sort"]
        for name in self._allowed_overrides:
            setattr(self, name, functools.partial(self._used_attribute, _name=name))

        if toclone:
            if custom_map is not None:
                self.extend(map(custom_map, toclone))
            elif plainstr:
                self.extend(map(str, toclone))
            elif strip_constraints:
                self.extend(map(strip_wildcard_constraints, toclone))
            else:
                self.extend(toclone)
            if isinstance(toclone, Namedlist):
                self._take_names(toclone._get_names())
        if fromdict:
            for key, item in fromdict.items():
                self.append(item)
                self._add_name(key)

    @staticmethod
    def _used_attribute(*args, _name, **kwargs):
        """
        Generic function that throws an `AttributeError`.

        Used as replacement for functions such as `index()` and `sort()`,
        which may be overridden by workflows, to signal to a user that
        these functions should not be used.
        """
        raise AttributeError(
            "{_name}() cannot be used; attribute name reserved"
            " for use in some existing workflows".format(_name=_name)
        )

    def _add_name(self, name):
        """
        Add a name to the last item.

        Arguments
        name -- a name
        """
        self._set_name(name, len(self) - 1)

    def _set_name(self, name, index, end=None):
        """
        Set the name of an item.

        Arguments
        name  -- a name
        index -- the item index
        """
        if name not in self._allowed_overrides and hasattr(self.__class__, name):
            raise AttributeError(
                "invalid name for input, output, wildcard, "
                "params or log: {name} is reserved for internal use".format(name=name)
            )

        self._names[name] = (index, end)
        if end is None:
            setattr(self, name, self[index])
        else:
            setattr(self, name, Namedlist(toclone=self[index:end]))

    def _get_names(self):
        """
        Get the defined names as (name, index) pairs.
        """
        for name, index in self._names.items():
            yield name, index

    def _take_names(self, names):
        """
        Take over the given names.

        Arguments
        names -- the given names as (name, index) pairs
        """
        for name, (i, j) in names:
            self._set_name(name, i, end=j)

    def items(self):
        for name in self._names:
            yield name, getattr(self, name)

    def _allitems(self):
        next = 0
        for name, index in sorted(
            self._names.items(),
            key=lambda item: (
                item[1][0],
                item[1][0] + 1 if item[1][1] is None else item[1][1],
            ),
        ):

            start, end = index
            if end is None:
                end = start + 1
            if start > next:
                for item in self[next:start]:
                    yield None, item
            yield name, getattr(self, name)
            next = end
        for item in self[next:]:
            yield None, item

    def _insert_items(self, index, items):
        self[index : index + 1] = items
        add = len(items) - 1
        for name, (i, j) in self._names.items():
            if i > index:
                self._names[name] = (i + add, None if j is None else j + add)
            elif i == index:
                self._set_name(name, i, end=i + len(items))

    def keys(self):
        return self._names.keys()

    def _plainstrings(self):
        return self.__class__.__call__(toclone=self, plainstr=True)

    def _stripped_constraints(self):
        return self.__class__.__call__(toclone=self, strip_constraints=True)

    def _clone(self):
        return self.__class__.__call__(toclone=self)

    def get(self, key, default_value=None):
        return self.__dict__.get(key, default_value)

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except TypeError:
            pass
        return getattr(self, key)

    def __hash__(self):
        return hash(tuple(self))

    def __str__(self):
        return " ".join(map(str, self))


class InputFiles(Namedlist):
    @property
    def size(self):
        return sum(f.size for f in self)

    @property
    def size_mb(self):
        return self.size / 1024 / 1024


class OutputFiles(Namedlist):
    pass


class Wildcards(Namedlist):
    pass


class Params(Namedlist):
    pass


class Resources(Namedlist):
    pass


class Log(Namedlist):
    pass


def _load_configfile(configpath_or_obj, filetype="Config"):
    "Tries to load a configfile first as JSON, then as YAML, into a dict."
    import yaml

    if isinstance(configpath_or_obj, str) or isinstance(configpath_or_obj, Path):
        obj = open(configpath_or_obj)
    else:
        obj = configpath_or_obj

    try:
        with obj as f:
            try:
                return json.load(f, object_pairs_hook=collections.OrderedDict)
            except ValueError:
                f.seek(0)  # try again
            try:
                # From https://stackoverflow.com/a/21912744/84349
                class OrderedLoader(yaml.Loader):
                    pass

                def construct_mapping(loader, node):
                    loader.flatten_mapping(node)
                    return collections.OrderedDict(loader.construct_pairs(node))

                OrderedLoader.add_constructor(
                    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping
                )
                return yaml.load(f, Loader=OrderedLoader)
            except yaml.YAMLError:
                raise WorkflowError(
                    "Config file is not valid JSON or YAML. "
                    "In case of YAML, make sure to not mix "
                    "whitespace and tab indentation.".format(filetype)
                )
    except FileNotFoundError:
        raise WorkflowError("{} file {} not found.".format(filetype, configpath))


def load_configfile(configpath):
    "Loads a JSON or YAML configfile as a dict, then checks that it's a dict."
    config = _load_configfile(configpath)
    if not isinstance(config, dict):
        raise WorkflowError(
            "Config file must be given as JSON or YAML " "with keys at top level."
        )
    return config


##### Wildcard pumping detection #####


class PeriodicityDetector:
    def __init__(self, min_repeat=20, max_repeat=100):
        """
        Args:
            max_repeat (int): The maximum length of the periodic substring.
            min_repeat (int): The minimum length of the periodic substring.
        """
        self.min_repeat = min_repeat
        self.regex = re.compile(
            "((?P<value>.+)(?P=value){{{min_repeat},{max_repeat}}})$".format(
                min_repeat=min_repeat - 1, max_repeat=max_repeat - 1
            )
        )

    def is_periodic(self, value):
        """Returns the periodic substring or None if not periodic."""
        # short-circuit: need at least min_repeat characters
        if len(value) < self.min_repeat:
            return None

        # short-circuit: need at least min_repeat same characters
        last_letter = value[-1]
        counter = collections.Counter(value)
        if counter[last_letter] < self.min_repeat:
            return None

        # short-circuit: need at least min_repeat same characters
        pos = 2
        while (
            value[-pos] != last_letter
        ):  # as long as last letter is not seen, repeat length is minimally pos
            if (
                len(value) < (pos * self.min_repeat)
                or counter[value[-pos]] < self.min_repeat
            ):
                return None
            pos += 1

        # now do the expensive regex
        m = self.regex.search(value)  # search for a periodic suffix.
        if m is not None:
            return m.group("value")
