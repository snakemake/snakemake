__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import ABC, abstractmethod
import asyncio
import collections
import copy
from dataclasses import dataclass, field
import datetime
import functools
import json
import os
import queue
import re
import shutil
import stat
import string
import subprocess as sp
import time
from contextlib import contextmanager
import string
import collections
import asyncio
from typing import Callable
from hashlib import sha256
from inspect import isfunction, ismethod
from itertools import chain, product
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Union

from snakemake_interface_common.utils import not_iterable
from snakemake_interface_storage_plugins.io import (
    WILDCARD_REGEX,
    IOCacheStorageInterface,
    Mtime,
    get_constant_prefix,
)

from snakemake.common import (
    ON_WINDOWS,
    async_run,
    get_function_params,
    get_input_function_aux_params,
    is_namedtuple_instance,
)
from snakemake.exceptions import (
    MissingOutputException,
    WildcardError,
    WorkflowError,
)
from snakemake.logging import logger


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


class AnnotatedStringInterface(ABC):
    @property
    @abstractmethod
    def flags(self) -> Dict[str, Any]:
        ...

    @abstractmethod
    def is_callable(self) -> bool:
        ...

    def is_flagged(self, flag: str) -> bool:
        return flag in self.flags and bool(self.flags[flag])


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
        if isinstance(path, _IOFile):
            parent = path.get_inventory_parent()
            return parent in self.has_inventory or super().__contains__(path)
        else:
            return super().__contains__(path)


class IOCache(IOCacheStorageInterface):
    def __init__(self, max_wait_time):
        self._mtime = dict()
        self._exists_local = ExistsDict(self)
        self._exists_in_storage = ExistsDict(self)
        self._size = dict()
        self.active = True
        self.remaining_wait_time = max_wait_time
        self.max_wait_time = max_wait_time

    @property
    def mtime(self):
        return self._mtime

    @property
    def exists_local(self):
        return self._exists_local

    @property
    def exists_in_storage(self):
        return self._exists_in_storage

    @property
    def size(self):
        return self._size

    async def mtime_inventory(self, jobs, n_workers=8):
        queue = asyncio.Queue()
        stop_item = object()

        async def worker(queue):
            while True:
                item = await queue.get()
                if item is stop_item:
                    queue.task_done()
                    return
                # Avoid superfluously checking mtime as the same file might be
                # added multiple times to the queue.
                if item not in self.mtime:
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
            for f in chain(job.input, job.output):
                if not f.is_storage and await f.exists():
                    queue.put_nowait(f)
            if (
                job.benchmark
                and not job.benchmark.is_storage
                and await job.benchmark.exists()
            ):
                queue.put_nowait(job.benchmark)

        # Send a stop item to each worker.
        for _ in range(n_workers):
            queue.put_nowait(stop_item)

        await asyncio.gather(*tasks)

    async def collect_mtime(self, path):
        return await path.mtime_uncached()

    def clear(self):
        self.mtime.clear()
        self.size.clear()
        self.exists_local.clear()
        self.exists_in_storage.clear()
        self.remaining_wait_time = self.max_wait_time

    def deactivate(self):
        self.clear()
        self.active = False


def IOFile(file, rule=None):
    assert rule is not None
    f = _IOFile(file)
    f.rule = rule
    return f


def iocache(func: Callable):
    @functools.wraps(func)
    async def wrapper(self, *args, **kwargs):
        # self: _IOFile
        if self.rule.workflow.iocache.active:
            cache = getattr(self.rule.workflow.iocache, func.__name__)
            if self in cache:
                return cache[self]
            v = await func(self, *args, **kwargs)
            cache[self] = v
            return v
        else:
            return await func(self, *args, **kwargs)

    return wrapper


class _IOFile(str, AnnotatedStringInterface):
    """
    A file that is either input or output of a rule.
    """

    __slots__ = ["_is_function", "_file", "rule", "_regex", "_wildcard_constraints"]

    def __new__(cls, file):
        is_annotated = isinstance(file, AnnotatedString)
        is_callable = (
            isfunction(file) or ismethod(file) or (is_annotated and bool(file.callable))
        )
        if isinstance(file, Path):
            file = str(file.as_posix())
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
        obj._wildcard_constraints = None

        if obj.is_storage:
            obj.storage_object._iofile = obj

        return obj

    def new_from(self, new_value):
        new = str.__new__(self.__class__, new_value)
        new._is_function = self._is_function
        new._file = self._file
        new.rule = self.rule
        if new.is_storage:
            new.storage_object._iofile = new
        return new

    async def inventory(self):
        """Starting from the given file, try to cache as much existence and
        modification date information of this and other files as possible.
        """
        cache = self.rule.workflow.iocache
        if cache.active:
            tasks = []
            if self.is_storage and self not in cache.exists_in_storage:
                # info not yet in inventory, let's discover as much as we can
                tasks.append(self.storage_object.inventory(cache))
            elif not ON_WINDOWS and self not in cache.exists_local:
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

        for i, path in enumerate(ancestors):
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

    def get_inventory_parent(self):
        """If eligible for inventory, get the parent of a given path.

        This code does not work on local Windows paths,
        but inventory is disabled on Windows.
        """
        if self.is_storage:
            return self.storage_object.get_inventory_parent()
        else:
            parent = os.path.dirname(self)
            if parent and parent != "..":
                return parent

    @contextmanager
    def open(self, mode="r", buffering=-1, encoding=None, errors=None, newline=None):
        """Open this file.

        This can (and should) be used in a `with`-statement.
        """
        f = open(self)
        try:
            yield f
        finally:
            f.close()

    def contains_wildcard(self):
        return contains_wildcard(self.file)

    @property
    def is_storage(self):
        return is_flagged(self._file, "storage_object")

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

    @property
    def should_keep_local(self):
        return self.storage_object.keep_local

    @property
    def should_not_be_retrieved_from_storage(self):
        return not self.storage_object.retrieve

    @property
    def storage_object(self):
        return get_flag_value(self._file, "storage_object")

    @property
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
        if _double_slash_regex.search(self._file) is not None and not self.is_storage:
            logger.warning(
                "File path {} contains double '{}'. "
                "This is likely unintended. {}".format(self._file, os.path.sep, hint)
            )

    async def exists(self):
        if self.is_storage:
            return await self.exists_in_storage()
        else:
            return await self.exists_local()

    def parents(self, omit=0):
        """Yield all parent paths, omitting the given number of ancestors."""
        for p in list(Path(self.file).parents)[::-1][omit:]:
            p = IOFile(str(p), rule=self.rule)
            p.clone_flags(self)
            yield p

    @iocache
    async def exists_local(self):
        return os.path.exists(self.file)

    @iocache
    async def exists_in_storage(self):
        if not self.is_storage:
            return False
        return await self.storage_object.managed_exists()

    async def protected(self):
        """Returns True if the file is protected. Always False for symlinks."""
        # symlinks are never regarded as protected
        return (
            await self.exists_local()
            and not os.access(self.file, os.W_OK)
            and not os.path.islink(self.file)
        )

    async def mtime(self):
        if self.rule.workflow.iocache.active:
            cache = self.rule.workflow.iocache
            if self in cache.mtime:
                mtime = cache.mtime[self]
                # if inventory is filled by storage plugin, mtime.local() will be None and
                # needs update
                if mtime.local() is None and await self.exists_local():
                    mtime_local = await self.mtime_uncached(skip_storage=True)
                    mtime._local_target = mtime_local._local_target
                    mtime._local = mtime_local._local
            else:
                cache.mtime[self] = mtime = await self.mtime_uncached()
            return mtime
        else:
            return await self.mtime_uncached()

    async def mtime_uncached(self, skip_storage: bool = False):
        """Obtain mtime.

        Usually, this will be one stat call only. For symlinks and directories
        it will be two, for symlinked directories it will be three,
        for storage files it will additionally query the storage
        location.
        """
        mtime_in_storage = (
            (await self.storage_object.managed_mtime())
            if self.is_storage and not skip_storage
            else None
        )

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
                return Mtime(local=mtime, storage=mtime_in_storage)

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
                    local=mtime, local_target=mtime_target, storage=mtime_in_storage
                )

        except FileNotFoundError:
            if self.is_storage:
                return Mtime(storage=mtime_in_storage)
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

    def is_fifo(self):
        """Return True if file is a FIFO according to the filesystem."""
        return stat.S_ISFIFO(os.stat(self).st_mode)

    @iocache
    async def size(self):
        if self.is_storage:
            try:
                return await self.storage_object.managed_size()
            except WorkflowError as e:
                try:
                    return await self.size_local()
                except IOError:
                    raise e
        else:
            return await self.size_local()

    async def size_local(self):
        # follow symlinks but throw error if invalid
        await self.check_broken_symlink()
        return os.path.getsize(self.file)

    async def is_checksum_eligible(self):
        return (
            await self.exists_local()
            and not os.path.isdir(self.file)
            and await self.size() < 100000
            and not self.is_fifo()
        )

    async def checksum(self, force=False):
        """Return checksum if file is small enough, else None.
        Returns None if file does not exist. If force is True,
        omit eligibility check."""
        if force or await self.is_checksum_eligible():  # less than 100000 bytes
            checksum = sha256()
            if await self.size() > 0:
                # only read if file is bigger than zero
                # otherwise the checksum is the same as taking hexdigest
                # from the empty sha256 as initialized above
                # This helps endless reading in case the input
                # is a named pipe or a socket or a symlink to a device like
                # /dev/random.
                with open(self.file, "rb") as f:
                    checksum.update(f.read())
            return checksum.hexdigest()
        else:
            return None

    async def is_same_checksum(self, other_checksum, force=False):
        checksum = await self.checksum(force=force)
        if checksum is None or other_checksum is None:
            # if no checksum available or files too large, not the same
            return False
        else:
            return checksum == other_checksum

    async def check_broken_symlink(self):
        """Raise WorkflowError if file is a broken symlink."""
        if not await self.exists_local():
            try:
                if os.lstat(self.file):
                    raise WorkflowError(
                        f"File {self.file} seems to be a broken symlink."
                    )
            except FileNotFoundError as e:
                # there is no broken symlink present, hence all fine
                return

    async def is_newer(self, time):
        """Returns true of the file (which is an input file) is newer than time, or if it is
        a symlink that points to a file newer than time."""
        if self.is_ancient:
            return False

        return (await self.mtime()).local_or_storage(follow_symlinks=True) > time

    async def retrieve_from_storage(self):
        if self.is_storage:
            if not self.should_not_be_retrieved_from_storage:

                async def is_newer_in_storage():
                    mtime = await self.mtime()
                    return mtime.local() < mtime.storage()

                if not await self.exists_local() or await is_newer_in_storage():
                    logger.info(f"Retrieving from storage: {self.storage_object.query}")
                    await self.storage_object.managed_retrieve()
                    logger.info("Finished retrieval.")
        else:
            raise WorkflowError(
                "The file to be retrieved does not seem to exist in the storage.",
                rule=self.rule,
            )

    async def store_in_storage(self):
        if self.is_storage:
            logger.info(f"Storing in storage: {self.storage_object.query}")
            await self.storage_object.managed_store()
            logger.info("Finished upload.")

    def prepare(self):
        dirpath = os.path.dirname(self.file)
        if len(dirpath) > 0:
            try:
                os.makedirs(dirpath, exist_ok=True)
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
        # iterate over content if output is a directory
        if os.path.isdir(self.file):
            # topdown=False ensures we chmod first the content, then the dir itself
            for dirpath, dirnames, filenames in os.walk(self.file, topdown=False):
                # no need to treat differently directories or files
                for content in dirnames + filenames:
                    lchmod(os.path.join(dirpath, content), mode)
        # protect explicit output itself
        lchmod(self.file, mode)

    async def remove(self, remove_non_empty_dir=False):
        if self.is_directory:
            await remove(self, remove_non_empty_dir=True)
        else:
            await remove(self, remove_non_empty_dir=remove_non_empty_dir)

    def touch(self, times=None):
        """times must be 2-tuple: (atime, mtime)"""
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

    def apply_wildcards(self, wildcards):
        f = self._file

        if self._is_function:
            f = self._file(Namedlist(fromdict=wildcards))

        # this bit ensures flags are transferred over to files after
        # wildcards are applied

        file_with_wildcards_applied = None

        if self.is_storage:
            query = apply_wildcards(self.storage_object.query, wildcards)
            storage_object = self.storage_object.__class__(
                query=query,
                provider=self.storage_object.provider,
                keep_local=self.storage_object.keep_local,
                retrieve=self.storage_object.retrieve,
            )

            validation_res = storage_object.is_valid_query()
            if not validation_res:
                raise WorkflowError(
                    str(validation_res),
                    rule=self.rule,
                )

            file_with_wildcards_applied = IOFile(
                storage_object.local_path(), rule=self.rule
            )
            file_with_wildcards_applied.clone_flags(self, skip_storage_object=True)
            file_with_wildcards_applied.flags["storage_object"] = storage_object
        else:
            file_with_wildcards_applied = IOFile(
                apply_wildcards(f, wildcards),
                rule=self.rule,
            )
            file_with_wildcards_applied.clone_flags(self)

        return file_with_wildcards_applied

    def get_wildcard_names(self):
        return get_wildcard_names(self.file)

    def regex(self):
        if self._regex is None:
            # compile a regular expression
            self._regex = re.compile(regex_from_filepattern(self.file))
        return self._regex

    def wildcard_constraints(self):
        if self._wildcard_constraints is None:
            self._wildcard_constraints = get_wildcard_constraints(self.file)
        return self._wildcard_constraints

    def constant_prefix(self):
        return get_constant_prefix(self.file)

    def constant_suffix(self):
        m = None
        for m in WILDCARD_REGEX.finditer(self.file):
            pass
        last_wildcard = m
        if last_wildcard:
            return self.file[last_wildcard.end() :]
        return self.file

    def match(self, target):
        return self.regex().match(target) or None

    def clone_flags(self, other, skip_storage_object=False):
        if isinstance(self._file, str):
            self._file = AnnotatedString(self._file)
        if isinstance(other._file, AnnotatedString) or isinstance(other._file, _IOFile):
            self._file.flags = getattr(other._file, "flags", {}).copy()
            if skip_storage_object and self.is_storage:
                del self._file.flags["storage_object"]
            else:
                assert "storage_object" not in self._file.flags, (
                    "bug: storage object cannot be safely cloned as the query might have "
                    "to change"
                )

    def clone_storage_object(self, other):
        if (
            isinstance(other._file, AnnotatedString)
            and "storage_object" in other._file.flags
        ):
            self._file.flags["storage_object"] = copy.copy(
                other._file.flags["storage_object"]
            )

    def set_flags(self, flags):
        if isinstance(self._file, str):
            self._file = AnnotatedString(self._file)
        self._file.flags = flags

    def __eq__(self, other):
        f = other._file if isinstance(other, _IOFile) else other
        return self._file == f

    def __hash__(self):
        return self._file.__hash__()


def pretty_print_iofile(iofile: _IOFile):
    if iofile.is_storage:
        return f"{iofile.storage_object.query} (storage)"
    else:
        return iofile._file


class AnnotatedString(str, AnnotatedStringInterface):
    def __init__(self, value):
        self._flags = dict()
        self.callable = value if is_callable(value) else None

    def new_from(self, new_value):
        new = str.__new__(self.__class__, new_value)
        new.flags = self.flags
        new.callable = self.callable
        return new

    def is_callable(self) -> bool:
        return self.callable is not None

    @property
    def flags(self) -> Dict[str, Any]:
        return self._flags

    @flags.setter
    def flags(self, value):
        self._flags = value


MaybeAnnotated = Union[AnnotatedStringInterface, str]


def is_flagged(value: MaybeAnnotated, flag: str) -> bool:
    if not isinstance(value, AnnotatedStringInterface):
        return False
    return value.is_flagged(flag)


def flag(value, flag_type, flag_value=True):
    if isinstance(value, AnnotatedStringInterface):
        value.flags[flag_type] = flag_value
        return value
    if not_iterable(value):
        if isinstance(value, Path):
            value = str(value.as_posix())
        value = AnnotatedString(value)
        value.flags[flag_type] = flag_value
        return value
    return [flag(v, flag_type, flag_value=flag_value) for v in value]


def is_callable(value: Any):
    return callable(value) or (
        isinstance(value, AnnotatedStringInterface) and value.is_callable()
    )


_double_slash_regex = (
    re.compile(r"([^:]//|^//)") if os.path.sep == "/" else re.compile(r"\\\\")
)


async def wait_for_files(
    files, latency_wait=3, wait_for_local=False, ignore_pipe_or_service=False
):
    """Wait for given files to be present in the filesystem."""
    files = list(files)

    async def get_missing():
        return [
            f
            for f in files
            if not (
                await f.exists_in_storage()
                if (
                    isinstance(f, _IOFile)
                    and f.is_storage
                    and (not wait_for_local or f.should_not_be_retrieved_from_storage)
                )
                else os.path.exists(f)
                if not (
                    (is_flagged(f, "pipe") or is_flagged(f, "service"))
                    and ignore_pipe_or_service
                )
                else True
            )
        ]

    missing = await get_missing()
    if missing:
        logger.info(f"Waiting at most {latency_wait} seconds for missing files.")
        for _ in range(latency_wait):
            missing = await get_missing()
            if not missing:
                return
            time.sleep(1)
        missing = "\n".join(await get_missing())
        raise IOError(
            f"Missing files after {latency_wait} seconds. This might be due to "
            "filesystem latency. If that is the case, consider to increase the "
            "wait time with --latency-wait:\n"
            f"{missing}"
        )


def get_wildcard_names(pattern):
    return set(match.group("name") for match in WILDCARD_REGEX.finditer(pattern))


def contains_wildcard(path):
    return WILDCARD_REGEX.search(str(path)) is not None


def contains_wildcard_constraints(pattern):
    return any(match.group("constraint") for match in WILDCARD_REGEX.finditer(pattern))


async def remove(file, remove_non_empty_dir=False):
    if file.is_storage and file.should_not_be_retrieved_from_storage:
        if await file.exists_in_storage():
            await file.storage_object.managed_remove()
    elif os.path.isdir(file) and not os.path.islink(file):
        if remove_non_empty_dir:
            shutil.rmtree(file)
        else:
            try:
                os.removedirs(file)
            except OSError as e:
                # skip non empty directories
                if e.errno == 39:
                    logger.info(f"Skipped removing non-empty directory {e.filename}")
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


def get_wildcard_constraints(pattern):
    constraints = {}
    for match in WILDCARD_REGEX.finditer(pattern):
        if match.group("constraint"):
            constraints[match.group("name")] = re.compile(match.group("constraint"))
    return constraints


def regex_from_filepattern(filepattern):
    f = []
    last = 0
    wildcards = set()
    for match in WILDCARD_REGEX.finditer(filepattern):
        f.append(re.escape(filepattern[last : match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "Constraint regex must be defined only in the first "
                    "occurrence of the wildcard in a string."
                )
            f.append(f"(?P={wildcard})")
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


def apply_wildcards(pattern, wildcards):
    def format_match(match):
        name = match.group("name")
        try:
            value = wildcards[name]
            return str(value)  # convert anything into a str
        except KeyError as ex:
            raise WildcardError(str(ex))

    return WILDCARD_REGEX.sub(format_match, pattern)


def is_callable(value):
    return (
        callable(value)
        or (isinstance(value, _IOFile) and value._is_function)
        or (isinstance(value, AnnotatedString) and value.callable is not None)
    )


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
    A flag to specify that output is a directory, rather than a file or named pipe.
    """
    if is_flagged(value, "pipe"):
        raise SyntaxError("Pipe and directory flags are mutually exclusive.")
    if is_flagged(value, "storage_object"):
        raise SyntaxError("Storage and directory flags are mutually exclusive.")
    return flag(value, "directory")


def temp(value):
    """
    A flag for an input or output file that shall be removed after usage.
    """
    if is_flagged(value, "protected"):
        raise SyntaxError("Protected and temporary flags are mutually exclusive.")
    if is_flagged(value, "storage_object"):
        raise SyntaxError("Storage and temporary flags are mutually exclusive.")
    return flag(value, "temp")


@dataclass
class QueueInfo:
    queue: queue.Queue
    finish_sentinel: Any
    last_checked: Optional[float] = None
    finished: bool = False
    items: List[Any] = field(default_factory=list, init=False)

    def consume(self, wildcards):
        assert (
            self.finished is False
        ), "bug: queue marked as finished but consume method called again"

        if wildcards:
            raise WorkflowError("from_queue() may not be used in rules with wildcards.")

        while True:
            try:
                item = self.queue.get_nowait()
            except queue.Empty:
                self.update_last_checked()
                return self.items
            self.queue.task_done()
            if item is self.finish_sentinel:
                logger.debug("finish sentinel found, stopping queue consumption")
                self.finished = True
                self.update_last_checked()
                return self.items
            self.items.append(item)

    def update_last_checked(self):
        self.last_checked = time.time()

    def __hash__(self):
        return hash(self.queue)


def from_queue(queue, finish_sentinel=None):
    if finish_sentinel is None:
        raise WorkflowError("Please provide a finish sentinel to from_queue.")

    queue_info = QueueInfo(queue, finish_sentinel)
    return flag(queue_info.consume, "from_queue", queue_info)


def pipe(value):
    if is_flagged(value, "protected"):
        raise SyntaxError("Pipes may not be protected.")
    if is_flagged(value, "storage_object"):
        raise SyntaxError("Pipes may not be in storage.")
    if ON_WINDOWS:
        logger.warning("Pipes are not yet supported on Windows.")
    return flag(value, "pipe", not ON_WINDOWS)


def service(value):
    if is_flagged(value, "protected"):
        raise SyntaxError("Service may not be protected.")
    if is_flagged(value, "storage_object"):
        raise SyntaxError("Service may not be storage files.")
    return flag(value, "service")


def temporary(value):
    """An alias for temp."""
    return temp(value)


def protected(value):
    """A flag for a file that shall be write-protected after creation."""
    if is_flagged(value, "temp"):
        raise SyntaxError("Protected and temporary flags are mutually exclusive.")
    if is_flagged(value, "storage_object"):
        raise SyntaxError("Storage and protected flags are mutually exclusive.")
    return flag(value, "protected")


def touch(value):
    return flag(value, "touch")


def ensure(value, non_empty=False, sha256=None):
    return flag(value, "ensure", {"non_empty": non_empty, "sha256": sha256})


def unpack(value):
    return flag(value, "unpack")


def repeat(value, n_repeat):
    """Flag benchmark records with the number of repeats."""
    return flag(value, "repeat", n_repeat)


def checkpoint_target(value):
    return flag(value, "checkpoint_target")


def sourcecache_entry(value, orig_path_or_uri):
    from snakemake.sourcecache import SourceFile

    assert not isinstance(
        orig_path_or_uri, SourceFile
    ), "bug: sourcecache_entry should receive a path or uri, not a SourceFile"
    return flag(value, "sourcecache_entry", orig_path_or_uri)


ReportObject = collections.namedtuple(
    "ReportObject",
    ["caption", "category", "subcategory", "labels", "patterns", "htmlindex"],
)


def report(
    value,
    caption=None,
    category=None,
    subcategory=None,
    labels=None,
    patterns=[],
    htmlindex=None,
):
    """Flag output file or directory as to be included into reports.

    In the case of a directory, files to include can be specified via a glob pattern (default: *).

    Arguments
    value -- File or directory.
    caption -- Path to a .rst file with a textual description of the result.
    category -- Name of the (optional) category in which the result should be displayed in the report.
    subcategory -- Name of the (optional) subcategory
    columns  -- Dict of strings (may contain wildcard expressions) that will be used as columns when displaying result tables
    patterns -- Wildcard patterns for selecting files if a directory is given (this is used as
               input for snakemake.io.glob_wildcards). Pattern shall not include the path to the
               directory itself.
    """
    return flag(
        value,
        "report",
        ReportObject(caption, category, subcategory, labels, patterns, htmlindex),
    )


def local(value):
    """Mark a file as a local file. This disables the application of a default storage
    provider.
    """
    if is_flagged(value, "storage_object"):
        raise SyntaxError("Storage and local flags are mutually exclusive.")
    return flag(value, "local")


def expand(*args, **wildcard_values):
    """
    Expand wildcards in given filepatterns.

    Arguments
    *args -- first arg: filepatterns as list or one single filepattern,
        second arg (optional): a function to combine wildcard values
        (itertools.product per default)
    **wildcard_values -- the wildcards as keyword arguments
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
    if "allow_missing" in wildcard_values and wildcard_values["allow_missing"] is True:

        class FormatDict(dict):
            def __missing__(self, key):
                return "{" + key + "}"

        format_dict = FormatDict
        # check that remove missing is not a wildcard in the filepatterns
        for filepattern in filepatterns:
            if "allow_missing" in re.findall(r"{([^}\.[!:]+)", filepattern):
                format_dict = dict
                break

    callables = {
        key: value
        for key, value in wildcard_values.items()
        if isinstance(value, Callable)
    }

    # remove unused wildcards to avoid duplicate filepatterns
    wildcard_values = {
        filepattern: {
            k: v
            for k, v in wildcard_values.items()
            if k in re.findall(r"{([^}\.[!:]+)", filepattern)
        }
        for filepattern in filepatterns
    }

    def do_expand(wildcard_values):
        def flatten(wildcard_values):
            for wildcard, values in wildcard_values.items():
                if (
                    isinstance(values, str)
                    or not isinstance(values, collections.abc.Iterable)
                    or is_namedtuple_instance(values)
                ):
                    values = [values]
                yield [(wildcard, value) for value in values]

        formatter = string.Formatter()
        try:
            return [
                formatter.vformat(filepattern, (), comb)
                for filepattern in filepatterns
                for comb in map(
                    format_dict, combinator(*flatten(wildcard_values[filepattern]))
                )
            ]
        except KeyError as e:
            raise WildcardError(f"No values given for wildcard {e}.")

    if callables:
        # defer expansion and return a function that does the expansion once it is called with
        # the usual arguments (wildcards, [input, ])
        def inner(wildcards, **aux_params):
            for wildcard, func in callables.items():
                func_aux_params = get_input_function_aux_params(func, aux_params)
                ret = func(wildcards, **func_aux_params)
                # store result for all filepatterns that need it
                for pattern_values in wildcard_values.values():
                    pattern_values[wildcard] = ret
            return do_expand(wildcard_values)

        return inner
    else:
        return do_expand(wildcard_values)


def multiext(prefix, *extensions):
    """Expand a given prefix with multiple extensions (e.g. .txt, .csv, _peaks.bed, ...)."""
    if any((r"/" in ext or r"\\" in ext) for ext in extensions):
        raise WorkflowError(
            r"Extensions for multiext may not contain path delimiters (/,\)."
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
    if is_flagged(pattern, "storage_object"):
        if files is not None:
            raise WorkflowError(
                "Error in glob_wildcards(): the files argument may not "
                "be combined with a storage object as pattern."
            )
        # for storage object patterns, we obtain the list of files from
        # the storage provider
        storage_object = pattern.flags["storage_object"]
        files = storage_object.list_candidate_matches()
        pattern = storage_object.query
    else:
        pattern = os.path.normpath(pattern)

    first_wildcard = re.search("{[^{]", pattern)
    dirname = (
        os.path.dirname(pattern[: first_wildcard.start()])
        if first_wildcard
        else os.path.dirname(pattern)
    )
    if not dirname:
        dirname = "."

    names = [match.group("name") for match in WILDCARD_REGEX.finditer(pattern)]
    Wildcards = collections.namedtuple("Wildcards", names)
    wildcards = Wildcards(*[list() for name in names])

    pattern = re.compile(regex_from_filepattern(pattern))

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
    pattern,
    wildcard_constraints: Dict[str, str],
    global_wildcard_constraints: Dict[str, str],
):
    """Update wildcard constraints

    Args:
      pattern (str): pattern on which to update constraints
      wildcard_constraints (dict): dictionary of wildcard:constraint key-value pairs
      global_wildcard_constraints (dict): dictionary of wildcard:constraint key-value pairs
    """

    def replace_constraint(match: re.Match):
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
            return f"{{{name},{newconstraint}}}"
        else:
            return match.group(0)

    examined_names: Set[str] = set()
    updated = WILDCARD_REGEX.sub(replace_constraint, pattern)

    # inherit flags
    if isinstance(pattern, AnnotatedString):
        updated = AnnotatedString(updated)
        updated.flags = dict(pattern.flags)
    return updated


def strip_wildcard_constraints(pattern):
    """Return a string that does not contain any wildcard constraints."""
    if is_callable(pattern):
        # do not apply on e.g. input functions
        return pattern

    def strip_constraint(match):
        return "{{{}}}".format(match.group("name"))

    return WILDCARD_REGEX.sub(strip_constraint, pattern)


class Namedlist(list):
    """
    A list that additionally provides functions to name items. Further,
    it is hashable, however, the hash does not consider the item names.
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

        if toclone is not None:
            if custom_map is not None:
                self.extend(map(custom_map, toclone))
            elif plainstr:
                self.extend(
                    # use original query if storage is not retrieved by snakemake
                    (str(x) if x.storage_object.retrieve else x.storage_object.query)
                    if isinstance(x, _IOFile) and x.storage_object is not None
                    else str(x)
                    for x in toclone
                )
            elif strip_constraints:
                self.extend(map(strip_wildcard_constraints, toclone))
            else:
                self.extend(toclone)
            if isinstance(toclone, Namedlist):
                self._take_names(toclone._get_names())
        if fromdict is not None:
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
        if isinstance(key, str):
            return getattr(self, key)
        else:
            return super().__getitem__(key)

    def __hash__(self):
        return hash(tuple(self))

    def __str__(self):
        return " ".join(map(str, self))


class InputFiles(Namedlist):
    @property
    def size(self):
        async def sizes():
            return [await f.size() for f in self]

        return sum(async_run(sizes()))

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
