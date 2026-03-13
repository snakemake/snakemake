import dataclasses
import json
from typing import Callable, Generic, Type, TypeVar, Any, overload

T = TypeVar("T")


def typed_to_dict(obj):
    """
    Convert NamedTuple/Dataclass/[objects with asdict function] to json.
    All the key must be valid varnames and all the values must be json serializable.
    """
    if dataclasses.is_dataclass(obj):
        if not isinstance(obj, type):
            return dataclasses.asdict(obj)
    elif isinstance(obj, tuple):
        if hasattr(obj, "_asdict"):  # NamedTuple
            return obj._asdict()  # type: ignore[reportAttributeAccessIssue]
    elif hasattr(obj, "asdict"):
        return obj.asdict()  # type: ignore[reportAttributeAccessIssue]
    raise NotImplementedError(f"Cannot convert {type(obj)} to dict")


class _TypedFile(str, Generic[T]):

    _loader: Callable[[str], T]
    _dumper: Callable[[T, str], None]
    def dump(self, obj: T, /):
        """Construct an instance of type_ and serialize it to the file as JSON."""
        self._dumper(obj, self)
    def load(self):
        """Deserialize the file and return an instance of type_."""
        return self._loader(self)  # type: ignore[arg-type]


class TypedFile(_TypedFile, Generic[T]):
    """
    A file path associated with a specific type for structured serialization.

    Behaves as a plain string (file path) in all contexts, but provides ``dump()`` and ``load()`` methods to serialize and deserialize instances of the associated type to and from disk as JSON.

    The associated type must be convertible to a dict (see :func:`typed_to_dict`).
    Type annotations serve as documentation only; no runtime type checking is performed.
    """
    type_: Type[T]
    def __new__(cls, value: str, type_: Type[T]):
        self = super().__new__(cls, value)
        self.type_ = type_
        return self
    def dump(self, *args, **kwargs):
        """Construct an instance of type_ and serialize it to the file as JSON."""
        obj = self.type_(*args, **kwargs)
        # obj1 = self.type_(**typed_to_dict(obj))  # is it needed to verify the type?
        self._dumper(obj, self)
    def load(self):
        """Deserialize the file and return an instance of type_."""
        return self.type_(**self._loader(self))  # type: ignore[arg-type]

@overload
def typed_factory(
        type_: None | str=None,
        *,
        loader: None | str | Callable[[str], T] =None,
        dumper: None | str | Callable[[T, str], None]=None
    ) -> Callable[[str], _TypedFile[T]]:...
@overload
def typed_factory(
        type_: Type[T],
        *,
        loader: None | str | Callable[[str], T] =None,
        dumper: None | str | Callable[[T, str], None]=None
    ) -> Callable[[str], TypedFile[T]]:...
def typed_factory(  # type: ignore[reportInconsistentOverloads]
        type_: None | str | Type[T]=None,
        *,
        loader: None | str | Callable[[str], T] =None,
        dumper: None | str | Callable[[T, str], None]=None
    ):
    """
    Returns a constructor for TypedFile bound to the given type.

    To determe how to load data from file, it could be:
    1. file type<str> => will be used to resolve the file format, and new a `T(**value)` instance.
    2. <None>         => file type will be inferred from actual file name, then follow 1.
    3. Callable[[str], type_] => ignore file type and replace default `.load(` method. Can be sth like `T.load` (classmethod called by `T.load(file)`) or `pd.read_csv`.

    To determe hwo to dump data to file, it could be:
    1. file type<str> => will be used to resolve the file format, then a `T(**value)` instance will be first generated.
                         For file types that should be converted to dict (e.g., json),
                         class `T` should be NamedTuple/Dataclass/[objects with asdict function],
                         and the `T(**value)` instance will be converted to dict before dumping.
    2. <None>         => file type will be inferred from actual file name (prior) or loader (only if it is a string), then follow 1.
    3. Callable[[T, str], None] => ignore file type and replace default `.dump(` method. Can be sth like `T.dump` (normally called by `T(sth).dump(file)`) or `pd.write_csv`.
    """
    typed_ = _TypedFile
    if type_ is None:
        assert loader or dumper
    elif isinstance(type_, str):
        loader = loader or type_
        dumper = dumper or type_
    else:
        typed_ = lambda file: TypedFile(file, type_)  # type: ignore[assignment]
    def _(file:str):
        typed_file: _TypedFile[T]|TypedFile[T] = typed_(file)
        loader_set = False
        if loader:
            loader_set = True
            if isinstance(loader, str):
                typed_file._loader, _dumper = resolve_file_format(loader)
            else:
                typed_file.load = lambda : loader(typed_file)  # type: ignore[method-assign]
        if dumper:
            if isinstance(dumper, str):
                if dumper == loader:
                    typed_file._dumper = _dumper
                else:
                    typed_file._dumper = resolve_file_format(dumper)[1]
            else:
                typed_file.dump = lambda obj: dumper(obj, typed_file)  # type: ignore[method-assign, misc]
        else:
            try:
                _loader, typed_file._dumper = resolve_file_format(file)
            except NotImplementedError:
                if isinstance(loader, str):
                    typed_file._dumper = _dumper
                else:
                    raise
            if not loader_set:
                typed_file._loader = _loader
        return typed_file
    return _


def resolve_file_format(file: str) -> tuple[Callable[[str], Any], Callable[[Any, str], None]]:
    suffix = file.rsplit(".", 1)[-1]
    open_: Callable[[str, str], Any]
    match suffix:
        case "gz":
            import gzip
            open_ = lambda f, mode: gzip.open(f, f"{mode}b")
        case "bz2":
            import bz2
            open_ = lambda f, mode: bz2.open(f, f"{mode}b")
        case "xz":
            import lzma
            open_ = lambda f, mode: lzma.open(f, f"{mode}b")
        case _:
            open_ = open
    if open_ is not open:
        suffix = file.rsplit(".", 2)[-2]
    match suffix:
        case "json":
            return (
                lambda f: json.load(open_(f, "r")),
                lambda obj, f: json.dump(typed_to_dict(obj), fp=open_(f, "w")),  # type: ignore
            )
        case "yaml" | "yml":
            import yaml
            return (
                lambda f: yaml.safe_load(open_(f, "r")),
                lambda obj, f: yaml.safe_dump(typed_to_dict(obj), stream=open_(f, "w")),  # type: ignore
            )
        case "toml":
            try:
                import toml
                return (
                    lambda f: toml.load(open_(f, "r")),  # type: ignore
                    lambda obj, f: toml.dump(typed_to_dict(obj), open_(f, "w")),  # type: ignore
                )
            except ImportError:
                def _(obj, f):
                    raise NotImplementedError("should use toml to dump toml file")
                import tomllib
                return (
                    lambda f: tomllib.load(open_(f, "rb")),  # type: ignore[arg-type]
                    _
                )
        case "pkl" | "pickle":
            import pickle
            return (
                lambda f: pickle.load(open_(f, "rb")),   # type: ignore[arg-type]
                lambda obj, f: pickle.dump(obj, open_(f, "wb")),  # type: ignore
            )
        case "npy":
            import numpy as np
            return (
                lambda f: np.load(f, allow_pickle=False),
                lambda obj, f: np.save(f, obj),
            )
        case "npz":
            import numpy as np
            return (
                lambda f: np.load(f, allow_pickle=False),
                lambda obj, f: np.savez_compressed(f, **typed_to_dict(obj)),
            )
        case "csv":
            import pandas as pd
            return (
                lambda f: pd.read_csv(f),
                lambda obj, f: obj.to_csv(f, index=False),  # type: ignore
            )
        case "tsv":
            import pandas as pd
            return (
                lambda f: pd.read_csv(f, sep="\t"),
                lambda obj, f: obj.to_csv(f, sep="\t", index=False),  # type: ignore
            )
        case "parquet":
            import pandas as pd
            return (
                lambda f: pd.read_parquet(f),
                lambda obj, f: obj.to_parquet(f, index=False),  # type: ignore
            )
        case "xlsx":
            import pandas as pd
            return (
                lambda f: pd.read_excel(f),
                lambda obj, f: obj.to_excel(f, index=False),  # type: ignore
            )
    raise NotImplementedError(f"Unsupported format[{suffix}] of file[{file}]")
