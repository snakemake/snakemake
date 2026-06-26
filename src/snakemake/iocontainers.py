__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

# Only modules from python standard library should be imported here (except for methods
# that are only called by Snakemake itself)!
# Further, the code here has to be kept compatible with Python 3.7.
# THe reason is that objects from this module are unpickled in scripts that might still
# run in older Python versions than Snakemake itself.

import os
import urllib.parse
from pathlib import Path

from typing import (
    Any,
    Callable,
    List,
    Dict,
    Union,
    Iterator,
    Tuple,
    TypeVar,
    Generic,
    Optional,
)

PathLike = Union[str, Path, os.PathLike]

FILE_HASH_PREFIX_LEN = 16


class ReportHref:
    def __init__(
        self,
        path: Union[str, Path],
        parent: Optional["ReportHref"] = None,
        url_args: Optional[Dict[str, str]] = None,
        anchor: Optional[str] = None,
    ):
        from snakemake.common import get_report_id

        self._parent = parent
        if parent is None:
            self._id = get_report_id(path)
        else:
            self._id = parent._id
        # ensure that path is a url compatible string
        self._path = path if isinstance(path, str) else str(path.as_posix())
        self._url_args = (
            {key: value for key, value in url_args.items()} if url_args else {}
        )
        self._anchor = anchor

    def child_path(self, path: Union[str, Path]):
        return self.__class__(path, parent=self)

    def url_args(self, **args: str):
        return self.__class__(path=self._path, parent=self._parent, url_args=args)

    def anchor(self, anchor: str):
        return self.__class__(
            path=self._path, parent=self._parent, url_args=self._url_args, anchor=anchor
        )

    def __str__(self) -> str:
        path = os.path.basename(self._path) if self._parent is None else self._path
        if self._url_args:

            def fmt_arg(key, value):
                return f"{key}={urllib.parse.quote(str(value))}"

            args = f"?{'&'.join(fmt_arg(key, value) for key, value in self._url_args.items())}"
        else:
            args = ""
        if self._anchor:
            anchor = f"#{urllib.parse.quote(self._anchor)}"
        else:
            anchor = ""
        return f"../{self._id[:FILE_HASH_PREFIX_LEN]}/{path}{args}{anchor}"


class AttributeGuard:
    def __init__(self, name):
        self.name = name

    def __call__(self, *args, **kwargs):
        """
        Generic function that throws an `AttributeError`.

        Used as replacement for functions such as `index()` and `sort()`,
        which may be overridden by workflows, to signal to a user that
        these functions should not be used.
        """
        raise AttributeError(
            f"{self.name}() cannot be used on snakemake input, output, resources etc.; "
            "instead it is a valid name for items on those objects. If you want e.g. to "
            "sort, convert to a plain list before or directly use sorted() on the "
            "object."
        )


# TODO: replace this with Self when Python 3.11 is the minimum supported version for
#   executing scripts
_TNamedList = TypeVar("_TNamedList")
"Type variable for self returning methods on Namedlist deriving classes"


class Namedlist(list, Generic[_TNamedList]):
    """
    A list that additionally provides functions to name items. Further,
    it is hashable, however, the hash does not consider the item names.
    """

    def __init__(
        self,
        toclone=None,
        fromdict: Optional[Dict[str, _TNamedList]] = None,
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
        from snakemake.io import _IOFile, strip_wildcard_constraints

        list.__init__(self)
        self._names: Dict[str, Tuple[int, int | None]] = dict()

        # white-list of attribute names that can be overridden in _set_name
        # default to throwing exception if called to prevent use as functions
        self._allowed_overrides = ["index", "sort"]
        for name in self._allowed_overrides:
            setattr(self, name, AttributeGuard(name))

        if toclone is not None:
            if custom_map is not None:
                self.extend(map(custom_map, toclone))
            elif plainstr:
                self.extend(
                    # use original query if storage is not retrieved by snakemake
                    (
                        (
                            str(x)
                            if x.storage_object.retrieve
                            else x.storage_object.query
                        )
                        if isinstance(x, _IOFile) and x.storage_object is not None
                        else str(x)
                    )
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

    def update(self, items: Dict):
        for key, value in items.items():
            if key in self._names:
                raise ValueError(f"Key {key} already exists in Namedlist")
            else:
                self.append(value)
                self._add_name(key)

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

    def items(self) -> Iterator[Tuple[str, _TNamedList]]:
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

    def _plainstrings(self: _TNamedList) -> _TNamedList:
        return self.__class__.__call__(toclone=self, plainstr=True)

    def _stripped_constraints(self: _TNamedList) -> _TNamedList:
        return self.__class__.__call__(toclone=self, strip_constraints=True)

    def _clone(self: _TNamedList) -> _TNamedList:
        return self.__class__.__call__(toclone=self)

    def get(self, key, default_value=None):
        value = self.__dict__.get(key, default_value)
        # handle internally guarded values like sort or index (see AttributeGuard)
        if isinstance(value, AttributeGuard):
            return default_value
        return value

    def __getitem__(self, key):
        if isinstance(key, str):
            return getattr(self, key)
        else:
            return super().__getitem__(key)

    def __hash__(self):  # type: ignore[override]
        return hash(tuple(self))

    def __str__(self):
        return " ".join(map(str, self))


class InputFiles(Namedlist):
    def _predicated_size_files(self, predicate: Callable) -> List[int]:
        from snakemake.common import async_run as async_run_fallback

        async def sizes() -> List[int]:
            import asyncio
            from snakemake.io import _IOFile

            async def get_size(f: _IOFile) -> Optional[int]:
                if await predicate(f):
                    return await f.size()
                return None

            sizes = await asyncio.gather(*map(get_size, self))
            return [res for res in sizes if res is not None]

        # the async_run method is expected to be provided through the eval context
        return globals().get("async_run", async_run_fallback)(sizes())

    @property
    def size_files(self):
        async def func_true(_) -> bool:
            return True

        return self._predicated_size_files(func_true)

    @property
    def size_tempfiles(self):
        from snakemake.io import is_flagged

        async def is_temp(iofile):
            return is_flagged(iofile, "temp")

        return self._predicated_size_files(is_temp)

    @property
    def size_files_kb(self):
        return [f / 1024 for f in self.size_files]

    @property
    def size_files_mb(self):
        return [f / 1024 for f in self.size_files_kb]

    @property
    def size_files_gb(self):
        return [f / 1024 for f in self.size_files_mb]

    @property
    def size(self):
        return sum(self.size_files)

    @property
    def size_kb(self):
        return sum(self.size_files_kb)

    @property
    def size_mb(self):
        return sum(self.size_files_mb)

    @property
    def temp_size_mb(self):
        return sum(self.size_tempfiles)

    @property
    def size_gb(self):
        return sum(self.size_files_gb)


class OutputFiles(Namedlist):
    pass


class Wildcards(Namedlist):
    pass


class Params(Namedlist):
    pass


class ResourceList(Namedlist):
    pass


class Log(Namedlist):
    pass


class Snakemake:
    def __init__(
        self,
        input_: InputFiles,
        output: OutputFiles,
        params: Params,
        wildcards: Wildcards,
        threads: int,
        resources: ResourceList,
        log: Log,
        config: Dict[str, Any],
        rulename: str,
        bench_iteration,
        scriptdir: Optional[PathLike] = None,
    ):
        # convert input and output to plain strings as some remote objects cannot
        # be pickled
        self.input = input_._plainstrings()
        self.output = output._plainstrings()
        self._safely_store_params(params)
        self.wildcards = wildcards
        self.threads = threads
        self.resources = resources
        self.log = log._plainstrings()
        self.config = config
        self.rule = rulename
        self.bench_iteration = bench_iteration
        self.scriptdir = scriptdir

    def report_href(self, path: Union[str, Path]) -> ReportHref:
        """Return an href to the given path in the report context, assuming that the
        path is given as it is given to the report marker in the workflow.

        The returned object can be extended to child paths using the `child_path(path)`
        method. This is useful if the referred item is a directory.
        """
        return ReportHref(path)

    def _log_shell_redirect(
        self,
        log: Optional[PathLike],
        stdout: bool = True,
        stderr: bool = True,
        append: bool = False,
    ) -> str:
        """
        Return a shell redirection string to be used in `shell()` calls

        This function allows scripts and wrappers support optional `log` files
        specified in the calling rule.  If no `log` was specified, then an
        empty string "" is returned, regardless of the values of `stdout`,
        `stderr`, and `append`.

        Parameters
        ---------

        stdout : bool
            Send stdout to log

        stderr : bool
            Send stderr to log

        append : bool
            Do not overwrite the log file. Useful for sending output of
            multiple commands to the same log. Note however that the log will
            not be truncated at the start.

        The following table describes the output:

        -------- -------- -------- ----- -------------
        stdout   stderr   append   log   return value
        -------- -------- -------- ----- ------------
        True     True     True     fn    >> fn 2>&1
        True     False    True     fn    >> fn
        False    True     True     fn    2>> fn
        True     True     False    fn    > fn 2>&1
        True     False    False    fn    > fn
        False    True     False    fn    2> fn
        any      any      any      None  ""
        -------- -------- -------- ----- -----------
        """
        if not log:
            return ""
        lookup = {
            (True, True, True): " >> {0} 2>&1",
            (True, False, True): " >> {0}",
            (False, True, True): " 2>> {0}",
            (True, True, False): " > {0} 2>&1",
            (True, False, False): " > {0}",
            (False, True, False): " 2> {0}",
        }
        return lookup[(stdout, stderr, append)].format(str(log))

    def log_fmt_shell(
        self, stdout: bool = True, stderr: bool = True, append: bool = False
    ) -> str:
        """
        Return a shell redirection string to be used in `shell()` calls

        This function allows scripts and wrappers to support optional `log` files
        specified in the calling rule.  If no `log` was specified, then an
        empty string "" is returned, regardless of the values of `stdout`,
        `stderr`, and `append`.

        Parameters
        ---------

        stdout : bool
            Send stdout to log

        stderr : bool
            Send stderr to log

        append : bool
            Do not overwrite the log file. Useful for sending an output of
            multiple commands to the same log. Note however that the log will
            not be truncated at the start.

        The following table describes the output:

        -------- -------- -------- ----- -------------
        stdout   stderr   append   log   return value
        -------- -------- -------- ----- ------------
        True     True     True     fn    >> fn 2>&1
        True     False    True     fn    >> fn
        False    True     True     fn    2>> fn
        True     True     False    fn    > fn 2>&1
        True     False    False    fn    > fn
        False    True     False    fn    2> fn
        any      any      any      None  ""
        -------- -------- -------- ----- -----------
        """
        return self._log_shell_redirect(str(self.log), stdout, stderr, append)

    @property
    def params(self):
        params = Params(toclone=list(self._params_store))
        try:
            for i, value in enumerate(params):
                param_type = self._params_types.get(i)
                if param_type is None:
                    # nothing to convert
                    continue
                if param_type.startswith("pd."):
                    import pandas as pd

                    if param_type == "pd.DataFrame":
                        params[i] = pd.DataFrame.from_dict(value)
                    elif param_type == "pd.Series":
                        params[i] = pd.Series(value)
                elif param_type.startswith("np."):
                    import numpy as np

                    if param_type == "np.ndarray":
                        params[i] = np.array(value)
                elif param_type.startswith("pl."):
                    import polars as pl

                    if param_type == "pl.LazyFrame":
                        params[i] = pl.from_dict(value).lazy()
                    elif param_type == "pl.DataFrame":
                        params[i] = pl.from_dict(value)
                    elif param_type == "pl.Series":
                        params[i] = pl.Series(**value)
        except ImportError as e:
            raise ImportError(
                "Failed to import required module for loading rule params. "
                "Make sure that the respective package (numpy, pandas, polars) "
                "is available in the software environment in which the "
                f"script/wrapper/notebook is executed: {e}"
            )

        params._take_names(self._params_store._get_names())
        return params

    def _safely_store_params(self, params: Params):
        try:
            import pandas as pd
        except ModuleNotFoundError:
            pd = None  # type: ignore[assignment]
        try:
            import numpy as np
        except ModuleNotFoundError:
            np = None  # type: ignore[assignment]
        try:
            import polars as pl
        except ModuleNotFoundError:
            pl = None  # type: ignore[assignment]

        self._params_store = Params(toclone=list(params))
        self._params_types = dict()
        for i, value in enumerate(params):
            if pd:
                if isinstance(value, pd.DataFrame):
                    self._params_store[i] = value.to_dict()
                    self._params_types[i] = "pd.DataFrame"
                elif isinstance(value, pd.Series):
                    self._params_store[i] = value.to_dict()
                    self._params_types[i] = "pd.Series"
            if np and isinstance(value, np.ndarray):
                self._params_store[i] = value.tolist()
                self._params_types[i] = "np.ndarray"
            if pl:
                if isinstance(value, pl.LazyFrame):
                    self._params_store[i] = value.collect().to_dict(as_series=False)
                    self._params_types[i] = "pl.LazyFrame"
                if isinstance(value, pl.DataFrame):
                    self._params_store[i] = value.to_dict(as_series=False)
                    self._params_types[i] = "pl.DataFrame"
                elif isinstance(value, pl.Series):
                    self._params_store[i] = {
                        "name": value.name,
                        "values": value.to_list(),
                    }
                    self._params_types[i] = "pl.Series"

        self._params_store._take_names(params._get_names())


# stub for the snakemake object, can be imported for type checking in scripts and wrappers
snakemake: Snakemake
