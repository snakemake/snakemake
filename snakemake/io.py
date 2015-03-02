# -*- coding: utf-8 -*-

import os
import re
import stat
import time
import json
from itertools import product, chain
from collections import Iterable, namedtuple
from snakemake.exceptions import MissingOutputException, WorkflowError, WildcardError
from snakemake.logging import logger

__author__ = "Johannes Köster"


def IOFile(file, rule=None):
    f = _IOFile(file)
    f.rule = rule
    return f


class _IOFile(str):
    """
    A file that is either input or output of a rule.
    """

    dynamic_fill = "__snakemake_dynamic__"

    def __new__(cls, file):
        obj = str.__new__(cls, file)
        obj._is_function = type(file).__name__ == "function"
        obj._file = file
        obj.rule = None
        obj._regex = None
        return obj

    @property
    def file(self):
        if not self._is_function:
            return self._file
        else:
            raise ValueError(
                "This IOFile is specified as a function and "
                "may not be used directly.")

    @property
    def exists(self):
        return os.path.exists(self.file)

    @property
    def protected(self):
        return self.exists and not os.access(self.file, os.W_OK)

    @property
    def mtime(self):
        return os.stat(self.file).st_mtime

    def is_newer(self, time):
        return self.mtime > time

    def prepare(self):
        path_until_wildcard = re.split(self.dynamic_fill, self.file)[0]
        dir = os.path.dirname(path_until_wildcard)
        if len(dir) > 0 and not os.path.exists(dir):
            try:
                os.makedirs(dir)
            except OSError as e:
                # ignore Errno 17 "File exists" (reason: multiprocessing)
                if e.errno != 17:
                    raise e

    def protect(self):
        mode = (os.stat(self.file).st_mode & ~stat.S_IWUSR &
            ~stat.S_IWGRP & ~stat.S_IWOTH)
        if os.path.isdir(self.file):
            for root, dirs, files in os.walk(self.file):
                for d in dirs:
                    os.chmod(os.path.join(self.file, d), mode)
                for f in files:
                    os.chmod(os.path.join(self.file, f), mode)
        else:
            os.chmod(self.file, mode)

    def remove(self):
        remove(self.file)

    def touch(self):
        try:
            os.utime(self.file, None)
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
            # create empty file
            with open(self.file, "w") as f:
                pass

    def apply_wildcards(
        self, wildcards, fill_missing=False,
        fail_dynamic=False):
        f = self._file
        if self._is_function:
            f = self._file(Namedlist(fromdict=wildcards))

        return IOFile(
            apply_wildcards(
                f, wildcards, fill_missing=fill_missing,
                fail_dynamic=fail_dynamic,
                dynamic_fill=self.dynamic_fill),
            rule=self.rule)

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

    def match(self, target):
        return self.regex().match(target) or None

    def format_dynamic(self):
        return self.replace(self.dynamic_fill, "{*}")

    def __eq__(self, other):
        f = other._file if isinstance(other, _IOFile) else other
        return self._file == f

    def __hash__(self):
        return self._file.__hash__()


_wildcard_regex = re.compile(
    "\{\s*(?P<name>\w+?)(\s*,\s*(?P<constraint>([^\{\}]+|\{\d+(,\d+)?\})*))?\s*\}")
#    "\{\s*(?P<name>\w+?)(\s*,\s*(?P<constraint>[^\}]*))?\s*\}")


def wait_for_files(files, latency_wait=3):
    """Wait for given files to be present in filesystem."""
    files = list(files)
    get_missing = lambda: [f for f in files if not os.path.exists(f)]
    missing = get_missing()
    if missing:
        logger.info("Waiting at most {} seconds for missing files.".format(latency_wait))
        for _ in range(latency_wait):
            if not get_missing():
                return
            time.sleep(1)
        raise IOError(
            "Missing files after {} seconds:\n{}".format(
                latency_wait, "\n".join(get_missing())))


def get_wildcard_names(pattern):
    return set(match.group('name') for match in
            _wildcard_regex.finditer(pattern))


def contains_wildcard(path):
    return _wildcard_regex.search(path) is not None


def remove(file):
    if os.path.exists(file):
        if os.path.isdir(file):
            try:
                os.removedirs(file)
            except OSError:
                # ignore non empty directories
                pass
        else:
            os.remove(file)


def regex(filepattern):
    f = []
    last = 0
    wildcards = set()
    for match in _wildcard_regex.finditer(filepattern):
        f.append(re.escape(filepattern[last:match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError("If multiple wildcards of the same name "
                "appear in a string, eventual constraints have to be defined "
                "at the first occurence and will be inherited by the others.")
            f.append("(?P={})".format(wildcard))
        else:
            wildcards.add(wildcard)
            f.append("(?P<{}>{})".format(
                wildcard,
                match.group("constraint")
                    if match.group("constraint") else ".+"))
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$") # ensure that the match spans the whole file
    return "".join(f)


def apply_wildcards(pattern, wildcards, fill_missing=False,
        fail_dynamic=False, dynamic_fill=None, keep_dynamic=False):

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
    return isinstance(value, str) or not isinstance(value, Iterable)


class AnnotatedString(str):
    def __init__(self, value):
        self.flags = dict()


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
        return flag in value.flags
    return False


def temp(value):
    """
    A flag for an input or output file that shall be removed after usage.
    """
    if is_flagged(value, "protected"):
        raise SyntaxError("Protected and temporary flags are mutually exclusive.")
    return flag(value, "temp")


def temporary(value):
    """ An alias for temp. """
    return temp(value)


def protected(value):
    """ A flag for a file that shall be write protected after creation. """
    if is_flagged(value, "temp"):
        raise SyntaxError("Protected and temporary flags are mutually exclusive.")
    return flag(value, "protected")


def dynamic(value):
    """
    A flag for a file that shall be dynamic, i.e. the multiplicity
    (and wildcard values) will be expanded after a certain
    rule has been run """
    annotated = flag(value, "dynamic")
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
        return [filepattern.format(**comb) for comb in map(dict, combinator(*flatten(wildcards))) for filepattern in filepatterns]
    except KeyError as e:
        raise WildcardError("No values given for wildcard {}.".format(e))


def glob_wildcards(pattern):
    """
    Glob the values of the wildcards by matching the given pattern to the filesystem.
    Returns a named tuple with a list of values for each wildcard.
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    dirname = os.path.dirname(pattern[:first_wildcard.start()]) if first_wildcard else os.path.dirname(pattern)
    if not dirname:
        dirname = "."
    
    names = [match.group('name')
        for match in _wildcard_regex.finditer(pattern)]
    Wildcards = namedtuple("Wildcards", names)
    wildcards = Wildcards(*[list() for name in names])

    pattern = re.compile(regex(pattern))
    for dirpath, dirnames, filenames in os.walk(dirname):
        for f in chain(filenames, dirnames):
            if dirpath != ".":
                f = os.path.join(dirpath, f)
            match = re.match(pattern, f)
            if match:
                for name, value in match.groupdict().items():
                    getattr(wildcards, name).append(value)
    return wildcards


# TODO rewrite Namedlist!
class Namedlist(list):
    """
    A list that additionally provides functions to name items. Further,
    it is hashable, however the hash does not consider the item names.
    """
    def __init__(self, toclone=None, fromdict=None, plainstr=False):
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
            self.extend(map(str, toclone) if plainstr else toclone)
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
        for name, index in sorted(
            self._names.items(),
            key=lambda item: item[1]):
            start, end = index
            if start > next:
                for item in self[next:start]:
                    yield None, item
            yield name, getattr(self, name)
            next = end if end is not None else start + 1
        for item in self[next:]:
            yield None, item

    def insert_items(self, index, items):
        self[index:index + 1] = items
        add = len(items) - 1
        for name, (i, j) in self._names.items():
            if i > index:
                self._names[name] = (i + add, j + add)
            elif i == index:
                self.set_name(name, i, end=i+len(items))

    def keys(self):
        return self._names

    def plainstrings(self):
        return self.__class__.__call__(toclone=self, plainstr=True)

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


def load_configfile(jsonpath):
    try:
        with open(jsonpath) as f:
            config = json.load(f)
    except FileNotFoundError:
        raise WorkflowError("Config file {} not found.".format(jsonpath))
    if not isinstance(config, dict):
        raise WorkflowError("Workflow config must be given as JSON with keys at top level.")
    return config


##### Wildcard pumping detection #####

class PeriodicityDetector:
    def __init__(self, min_repeat=50, max_repeat=100):
        """
        Args:
            max_len (int): The maximum length of the periodic substring.
        """
        self.regex = re.compile("((?P<value>.+)(?P=value){{{min_repeat},{max_repeat}}})$".format(min_repeat=min_repeat - 1, max_repeat=max_repeat - 1))

    def is_periodic(self, value):
        """Returns the periodic substring or None if not periodic."""
        m = self.regex.search(value) # search for a periodic suffix.
        if m is not None:
            return m.group("value")
