__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import collections
import git
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
from collections import Iterable, namedtuple
from snakemake.exceptions import MissingOutputException, WorkflowError, WildcardError, RemoteFileException
from snakemake.logging import logger
from inspect import isfunction, ismethod

from snakemake.common import DYNAMIC_FILL


def lstat(f):
    return os.stat(f,
                   follow_symlinks=os.stat not in os.supports_follow_symlinks)


def lutime(f, times):
    #In some cases, we have a platform where os.supports_follow_symlink includes stat()
    #but not utime().  This leads to an anomaly.  In any case we never want to touch the
    #target of a link.
    if os.utime in os.supports_follow_symlinks:
        #...utime is well behaved
        os.utime(f, times, follow_symlinks=False)
    elif not os.path.islink(f):
        #...symlinks not an issue here
        os.utime(f, times)
    else:
        try:
            # try the system command
            if times:
                fmt_time = lambda sec: datetime.fromtimestamp(sec).strftime("%Y%m%d%H%M.%S")
                atime, mtime = times
                sp.check_call(["touch", "-h", f, "-a", "-t", fmt_time(atime)])
                sp.check_call(["touch", "-h", f, "-m", "-t", fmt_time(mtime)])
            else:
                sp.check_call(["touch", "-h", f])
        except sp.CalledProcessError:
            pass
        #...problem system.  Do nothing.
        logger.warning("Unable to set utime on symlink {}. Your Python build does not support it.".format(f))
        return None


def lchmod(f, mode):
    os.chmod(f,
             mode,
             follow_symlinks=os.chmod not in os.supports_follow_symlinks)


class IOCache:
    def __init__(self):
        self.mtime = dict()
        self.exists_local = dict()
        self.exists_remote = dict()
        self.size = dict()
        self.active = True

    def clear(self):
        self.mtime.clear()
        self.size.clear()
        self.exists_local.clear()
        self.exists_remote.clear()

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

    __slots__ = ["_is_function", "_file", "rule", "_regex"]

    def __new__(cls, file):
        obj = str.__new__(cls, file)
        obj._is_function = isfunction(file) or ismethod(file)
        obj._is_function = obj._is_function or (
            isinstance(file, AnnotatedString) and bool(file.callable))
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
                    return getattr(self.remote_object, func.__name__)(*args, **
                                                                      kwargs)
            return func(self, *args, **kwargs)

        return wrapper

    @property
    def is_remote(self):
        return is_flagged(self._file, "remote_object")

    @property
    def is_ancient(self):
        return is_flagged(self._file, "ancient")

    @property
    def is_directory(self):
        return is_flagged(self._file, "directory")

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
            raise ValueError("This IOFile is specified as a function and "
                             "may not be used directly.")

    def check(self):
        hint = (
            "It can also lead to inconsistent results of the file-matching "
            "approach used by Snakemake."
        )
        if self._file.startswith("./"):
            logger.warning("Relative file path '{}' starts with './'. This is redundant "
                           "and strongly discouraged. {} You can simply omit the './' "
                           "for relative file paths.".format(self._file, hint))
        if self._file.startswith(" "):
            logger.warning("File path '{}' starts with whitespace. "
                "This is likely unintended. {}".format(self._file, hint))
        if self._file.endswith(" "):
            logger.warning("File path '{}' ends with whitespace. "
                "This is likely unintended. {}".format(self._file, hint))
        if "\n" in self._file:
            logger.warning("File path '{}' contains line break. "
                "This is likely unintended. {}".format(self._file, hint))
        if _double_slash_regex.search(self._file) is not None:
            logger.warning("File path {} contains double '{}'. "
                "This is likely unintended. {}".format(
                    self._file, os.path.sep, hint))

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
        if self.rule.workflow.iocache.active:
            # The idea is to first check existence of parent directories and
            # cache the results.
            # We omit the last ancestor, because this is always "." or "/" or a
            # drive letter.
            for p in self.parents(omit=1):
                try:
                    if not p.exists_local:
                        return False
                except:
                    # In case of an error, we continue, because it can be that
                    # we simply don't have the permissions to access a parent
                    # directory.
                    continue
        return os.path.exists(self.file)

    @property
    @iocache
    def exists_remote(self):
        if not self.is_remote:
            return False
        if (self.rule.workflow.iocache.active and
            self.remote_object.provider.allows_directories):
            # The idea is to first check existence of parent directories and
            # cache the results.
            # We omit the last 2 ancestors, because these are "." and the host
            # name of the remote location.
            for p in self.parents(omit=2):
                try:
                    if not p.exists_remote:
                        return False
                except:
                    # In case of an error, we continue, because it can be that
                    # we simply don't have the permissions to access a parent
                    # directory in the remote.
                    continue
        return self.remote_object.exists()

    @property
    def protected(self):
        return self.exists_local and not os.access(self.file, os.W_OK)

    @property
    @iocache
    @_refer_to_remote
    def mtime(self):
        return self.mtime_local

    @property
    def mtime_local(self):
        # do not follow symlinks for modification time
        if os.path.isdir(self.file) and os.path.exists(os.path.join(self.file, ".snakemake_timestamp")):
            return lstat(os.path.join(self.file, ".snakemake_timestamp")).st_mtime
        else:
            return lstat(self.file).st_mtime

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
        if not self.exists_local and lstat(self.file):
            raise WorkflowError("File {} seems to be a broken symlink.".format(
                self.file))

    @_refer_to_remote
    def is_newer(self, time):
        """ Returns true of the file is newer than time, or if it is
            a symlink that points to a file newer than time. """
        if self.is_ancient:
            return False
        elif self.is_remote:
            #If file is remote but provider does not override the implementation this
            #is the best we can do.
            return self.mtime > time
        else:
            if os.path.isdir(self.file) and os.path.exists(os.path.join(self.file, ".snakemake_timestamp")):
                st_mtime_file = os.path.join(self.file, ".snakemake_timestamp")
            else:
                st_mtime_file = self.file
            try:
                return os.stat(st_mtime_file, follow_symlinks=True).st_mtime > time or self.mtime > time
            except FileNotFoundError:
                raise WorkflowError("File {} not found although it existed before. Is there another active process that might have deleted it?")

    def download_from_remote(self):
        if self.is_remote and self.remote_object.exists():
            if not self.should_stay_on_remote:
                logger.info("Downloading from remote: {}".format(self.file))
                self.remote_object.download()
                logger.info("Finished download.")
        else:
            raise RemoteFileException(
                "The file to be downloaded does not seem to exist remotely.")

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
                os.makedirs(dir, exist_ok = True)
            except OSError as e:
                # ignore Errno 17 "File exists" (reason: multiprocessing)
                if e.errno != 17:
                    raise e

        if is_flagged(self._file, "pipe"):
            os.mkfifo(self._file)

    def protect(self):
        mode = (lstat(self.file).st_mode & ~stat.S_IWUSR & ~stat.S_IWGRP
                & ~stat.S_IWOTH)
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
                    with open(file, "w") as f:
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
                    snakefile=self.rule.snakefile)
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
            file = os.path.join(self.file, ".snakemake_timestamp") if self.is_directory else self.file
            with open(file, "w") as f:
                pass

    def apply_wildcards(self,
                        wildcards,
                        fill_missing=False,
                        fail_dynamic=False):
        f = self._file
        if self._is_function:
            f = self._file(Namedlist(fromdict=wildcards))

        # this bit ensures flags are transferred over to files after
        # wildcards are applied

        file_with_wildcards_applied = IOFile(
            apply_wildcards(f,
                            wildcards,
                            fill_missing=fill_missing,
                            fail_dynamic=fail_dynamic,
                            dynamic_fill=DYNAMIC_FILL),
            rule=self.rule)

        file_with_wildcards_applied.clone_flags(self)

        return file_with_wildcards_applied

    def get_wildcard_names(self):
        return get_wildcard_names(self.file)

    def contains_wildcard(self):
        return contains_wildcard(self.file)

    def regex(self):
        if self._regex is None:
            # compile a regular expression
            self._regex = re.compile(regex(self.file))
        return self._regex

    def constant_prefix(self):
        first_wildcard = _wildcard_regex.search(self.file)
        if first_wildcard:
            return self.file[:first_wildcard.start()]
        return self.file

    def constant_suffix(self):
        m = None
        for m in _wildcard_regex.finditer(self.file):
            pass
        last_wildcard = m
        if last_wildcard:
            return self.file[last_wildcard.end():]
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
                self._file.flags['remote_object'] = copy.copy(
                    self._file.flags['remote_object'])
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


_double_slash_regex = (re.compile(r"([^:]//|^//)")
                       if os.path.sep == "/"
                       else re.compile(r"\\\\"))


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
    """, re.VERBOSE)


def wait_for_files(files, latency_wait=3, force_stay_on_remote=False):
    """Wait for given files to be present in filesystem."""
    files = list(files)
    def get_missing():
        return [
            f for f in files
            if not (f.exists_remote
                    if (isinstance(f, _IOFile) and
                       f.is_remote and
                       (force_stay_on_remote or f.should_stay_on_remote))
                    else os.path.exists(f))]

    missing = get_missing()
    if missing:
        logger.info("Waiting at most {} seconds for missing files.".format(
            latency_wait))
        for _ in range(latency_wait):
            if not get_missing():
                return
            time.sleep(1)
        raise IOError("Missing files after {} seconds:\n{}".format(
            latency_wait, "\n".join(get_missing())))


def get_wildcard_names(pattern):
    return set(match.group('name')
               for match in _wildcard_regex.finditer(pattern))


def contains_wildcard(path):
    return _wildcard_regex.search(path) is not None


def contains_wildcard_constraints(pattern):
    return any(match.group('constraint') for match in _wildcard_regex.finditer(pattern))


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
                    logger.info("Skipped removing non-empty directory {}".format(e.filename))
                else:
                    logger.warning(str(e))
    #Remember that dangling symlinks fail the os.path.exists() test, but
    #we definitely still want to zap them. try/except is the safest way.
    #Also, we don't want to remove the null device if it is an output.
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
        f.append(re.escape(filepattern[last:match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "Constraint regex must be defined only in the first "
                    "occurence of the wildcard in a string.")
            f.append("(?P={})".format(wildcard))
        else:
            wildcards.add(wildcard)
            f.append("(?P<{}>{})".format(wildcard, match.group("constraint") if
                                         match.group("constraint") else ".+"))
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$")  # ensure that the match spans the whole file
    return "".join(f)


def apply_wildcards(pattern,
                    wildcards,
                    fill_missing=False,
                    fail_dynamic=False,
                    dynamic_fill=None,
                    keep_dynamic=False):
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

    return re.sub(_wildcard_regex, format_match, pattern)


def not_iterable(value):
    return isinstance(value, str) or isinstance(value, dict) or not isinstance(
        value, Iterable)


def is_callable(value):
    return (callable(value) or
            (isinstance(value, _IOFile) and value._is_function))


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
        raise SyntaxError(
            "Protected and temporary flags are mutually exclusive.")
    if is_flagged(value, "remote"):
        raise SyntaxError("Remote and temporary flags are mutually exclusive.")
    return flag(value, "temp")


def pipe(value):
    if is_flagged(value, "protected"):
        raise SyntaxError("Pipes may not be protected.")
    if is_flagged(value, "remote"):
        raise SyntaxError("Pipes may not be remote files.")
    return flag(value, "pipe")


def temporary(value):
    """ An alias for temp. """
    return temp(value)


def protected(value):
    """ A flag for a file that shall be write protected after creation. """
    if is_flagged(value, "temp"):
        raise SyntaxError(
            "Protected and temporary flags are mutually exclusive.")
    if is_flagged(value, "remote"):
        raise SyntaxError("Remote and protected flags are mutually exclusive.")
    return flag(value, "protected")


def dynamic(value):
    """
    A flag for a file that shall be dynamic, i.e. the multiplicity
    (and wildcard values) will be expanded after a certain
    rule has been run """
    annotated = flag(value, "dynamic", True)
    tocheck = [annotated] if not_iterable(annotated) else annotated
    for file in tocheck:
        matches = list(_wildcard_regex.finditer(file))
        #if len(matches) != 1:
        #    raise SyntaxError("Dynamic files need exactly one wildcard.")
        for match in matches:
            if match.group("constraint"):
                raise SyntaxError(
                    "The wildcards in dynamic files cannot be constrained.")
    return annotated


def touch(value):
    return flag(value, "touch")


def unpack(value):
    return flag(value, "unpack")


def repeat(value, n_repeat):
    """Flag benchmark records with the number of repeats."""
    return flag(value, "repeat", n_repeat)


ReportObject = namedtuple("ReportObject", ["caption", "category"])


def report(value, caption=None, category=None):
    """Flag output file as to be included into reports."""
    return flag(value, "report", ReportObject(caption, category))


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
        with their values as lists
    """
    filepatterns = args[0]
    if len(args) == 1:
        combinator = product
    elif len(args) == 2:
        combinator = args[1]
    if isinstance(filepatterns, str):
        filepatterns = [filepatterns]

    def flatten(wildcards):
        for wildcard, values in wildcards.items():
            if isinstance(values, str) or not isinstance(values, Iterable):
                values = [values]
            yield [(wildcard, value) for value in values]

    try:
        return [filepattern.format(**comb)
                for comb in map(dict, combinator(*flatten(wildcards)))
                for filepattern in filepatterns]
    except KeyError as e:
        raise WildcardError("No values given for wildcard {}.".format(e))


def limit(pattern, **wildcards):
    """
    Limit wildcards to the given values.

    Arguments:
    **wildcards -- the wildcards as keyword arguments
                   with their values as lists
    """
    return pattern.format(**{
        wildcard: "{{{},{}}}".format(wildcard, "|".join(values))
        for wildcard, values in wildcards.items()
    })


def glob_wildcards(pattern, files=None):
    """
    Glob the values of the wildcards by matching the given pattern to the filesystem.
    Returns a named tuple with a list of values for each wildcard.
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    dirname = os.path.dirname(pattern[:first_wildcard.start(
    )]) if first_wildcard else os.path.dirname(pattern)
    if not dirname:
        dirname = "."

    names = [match.group('name')
             for match in _wildcard_regex.finditer(pattern)]
    Wildcards = namedtuple("Wildcards", names)
    wildcards = Wildcards(*[list() for name in names])

    pattern = re.compile(regex(pattern))

    if files is None:
        files = (os.path.normpath(os.path.join(dirpath, f))
                 for dirpath, dirnames, filenames in os.walk(dirname)
                 for f in chain(filenames, dirnames))

    for f in files:
        match = re.match(pattern, f)
        if match:
            for name, value in match.groupdict().items():
                getattr(wildcards, name).append(value)
    return wildcards


def update_wildcard_constraints(pattern,
                                wildcard_constraints,
                                global_wildcard_constraints):
    """Update wildcard constraints

    Args:
      pattern (str): pattern on which to update constraints
      wildcard_constraints (dict): dictionary of wildcard:constraint key-value pairs
      global_wildcard_constraints (dict): dictionary of wildcard:constraint key-value pairs
    """
    def replace_constraint(match):
        name = match.group("name")
        constraint = match.group("constraint")
        newconstraint = wildcard_constraints.get(name, global_wildcard_constraints.get(name))
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
    updated = re.sub(_wildcard_regex, replace_constraint, pattern)

    # inherit flags
    if isinstance(pattern, AnnotatedString):
        updated = AnnotatedString(updated)
        updated.flags = dict(pattern.flags)
    return updated

def split_git_path(path):
    file_sub = re.sub("^git\+file:/+",'/',path)
    (file_path, version) = file_sub.split("@")
    file_path = os.path.realpath(file_path)
    root_path = get_git_root(file_path)
    if file_path.startswith(root_path):
        file_path = file_path[len(root_path):].lstrip("/")
    return (root_path, file_path, version)


def get_git_root(path):
    """
        Args:
            path: (str) Path a to a directory/file that is located inside the repo
        Returns:
            path to root folder for git repo
    """
    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        return git_repo.git.rev_parse("--show-toplevel")
    except git.exc.NoSuchPathError as e:
        tail,head = os.path.split(path)
        return get_git_root_parent_directory(tail,path)

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
    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        return git_repo.git.rev_parse("--show-toplevel")
    except git.exc.NoSuchPathError as e:
        tail,head = os.path.split(path)
        if tail is None:
            raise WorkflowError("Neither provided git path ({}) ".format(input_path) +
                                "or parent directories contain a valid git repo.")
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
    if git_file.startswith("git+file:"):
        (root_path, file_path, version) = split_git_path(git_file)
        return git.Repo(root_path).git.show('{}:{}'.format(version, file_path))
    else:
        raise WorkflowError("Provided git path ({}) doesn't meet the "
                        "expected format:".format(git_file) +
                        ", expected format is "
                        "git+file://PATH_TO_REPO/PATH_TO_FILE_INSIDE_REPO@VERSION")
    return None

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

    def __init__(self, toclone=None, fromdict=None,
                 plainstr=False, strip_constraints=False):
        """
        Create the object.

        Arguments
        toclone  -- another Namedlist that shall be cloned
        fromdict -- a dict that shall be converted to a
            Namedlist (keys become names)
        """
        list.__init__(self)
        self._names = dict()

        if toclone:
            if plainstr:
                self.extend(map(str, toclone))
            elif strip_constraints:
                self.extend(map(strip_wildcard_constraints, toclone))
            else:
                self.extend(toclone)
            if isinstance(toclone, Namedlist):
                self.take_names(toclone.get_names())
        if fromdict:
            for key, item in fromdict.items():
                self.append(item)
                self.add_name(key)

    def add_name(self, name):
        """
        Add a name to the last item.

        Arguments
        name -- a name
        """
        self.set_name(name, len(self) - 1)

    def set_name(self, name, index, end=None):
        """
        Set the name of an item.

        Arguments
        name  -- a name
        index -- the item index
        """
        self._names[name] = (index, end)
        if end is None:
            setattr(self, name, self[index])
        else:
            setattr(self, name, Namedlist(toclone=self[index:end]))

    def get_names(self):
        """
        Get the defined names as (name, index) pairs.
        """
        for name, index in self._names.items():
            yield name, index

    def take_names(self, names):
        """
        Take over the given names.

        Arguments
        names -- the given names as (name, index) pairs
        """
        for name, (i, j) in names:
            self.set_name(name, i, end=j)

    def items(self):
        for name in self._names:
            yield name, getattr(self, name)

    def allitems(self):
        next = 0
        for name, index in sorted(self._names.items(),
                key=lambda item: (item[1][0], item[1][0] + 1 if item[1][1] is None else item[1][1])):

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

    def insert_items(self, index, items):
        self[index:index + 1] = items
        add = len(items) - 1
        for name, (i, j) in self._names.items():
            if i > index:
                self._names[name] = (i + add, None if j is None else j + add)
            elif i == index:
                self.set_name(name, i, end=i + len(items))

    def keys(self):
        return self._names

    def plainstrings(self):
        return self.__class__.__call__(toclone=self, plainstr=True)

    def stripped_constraints(self):
        return self.__class__.__call__(toclone=self, strip_constraints=True)

    def clone(self):
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
    pass


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


def _load_configfile(configpath, filetype="Config"):
    "Tries to load a configfile first as JSON, then as YAML, into a dict."
    try:
        with open(configpath) as f:
            try:
                return json.load(f, object_pairs_hook=collections.OrderedDict)
            except ValueError:
                f.seek(0)  # try again
            try:
                import yaml
            except ImportError:
                raise WorkflowError("{} file is not valid JSON and PyYAML "
                                    "has not been installed. Please install "
                                    "PyYAML to use YAML config files.".format(
                                    filetype))
            try:
                # From http://stackoverflow.com/a/21912744/84349
                class OrderedLoader(yaml.Loader):
                    pass
                def construct_mapping(loader, node):
                    loader.flatten_mapping(node)
                    return collections.OrderedDict(
                        loader.construct_pairs(node))
                OrderedLoader.add_constructor(
                    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                    construct_mapping)
                return yaml.load(f, OrderedLoader)
            except yaml.YAMLError:
                raise WorkflowError("Config file is not valid JSON or YAML. "
                                    "In case of YAML, make sure to not mix "
                                    "whitespace and tab indentation.".format(
                                    filetype))
    except FileNotFoundError:
        raise WorkflowError("{} file {} not found.".format(filetype, configpath))


def load_configfile(configpath):
    "Loads a JSON or YAML configfile as a dict, then checks that it's a dict."
    config = _load_configfile(configpath)
    if not isinstance(config, dict):
        raise WorkflowError("Config file must be given as JSON or YAML "
                            "with keys at top level.")
    return config

##### Wildcard pumping detection #####


class PeriodicityDetector:
    def __init__(self, min_repeat=20, max_repeat=100):
        """
        Args:
            max_repeat (int): The maximum length of the periodic substring.
            min_repeat (int): The minimum length of the periodic substring.
        """
        self.regex = re.compile(
            "((?P<value>.+)(?P=value){{{min_repeat},{max_repeat}}})$".format(
                min_repeat=min_repeat - 1,
                max_repeat=max_repeat - 1))

    def is_periodic(self, value):
        """Returns the periodic substring or None if not periodic."""
        m = self.regex.search(value)  # search for a periodic suffix.
        if m is not None:
            return m.group("value")
