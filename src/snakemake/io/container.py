from __future__ import annotations

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from typing import (
    Callable,
    List,
    TypeVar,
    Generic,
)


# TODO: replace this with Self when Python 3.11 is the minimum supported version for
#   executing scripts
_TNamedList = TypeVar("_TNamedList")
"Type variable for self returning methods on Namedlist deriving classes"

_TNamedKeys = TypeVar("_TNamedKeys")
"Type variable for self returning methods on Namedlist deriving classes"

class Namedlist(list, Generic[_TNamedKeys, _TNamedList]):
    """
    A list that additionally provides functions to name items. Further,
    it is hashable, however, the hash does not consider the item names.
    """

    def __init__(
        self,
        toclone=None,
        fromdict: Optional[Dict[_TNamedKeys, _TNamedList]] = None,
        plainstr=False,
        strip_constraints=False,
        custom_map=None,
    ):
        from snakemake.io import _IOFile, AttributeGuard
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

    def items(self) -> Iterator[Tuple[_TNamedKeys, _TNamedList]]:
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
        from snakemake.io import AttributeGuard
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

    def __hash__(self):
        return hash(tuple(self))

    def __str__(self):
        return " ".join(map(str, self))


class InputFiles(Namedlist):
    def _predicated_size_files(self, predicate: Callable) -> List[int]:
        async def sizes() -> List[int]:
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
