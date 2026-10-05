from typing import Callable, Generic, Type, TypeVar, Any, overload, Tuple
import sys

T = TypeVar("T")


def disabled_func(disabled: str, supplied: str):
    def _(*args, **kwargs):
        raise ValueError(
            f"{disabled} is disabled when {supplied} is supplied without {disabled}. "
        )

    return _


class _TypedFile(str, Generic[T]):

    _loader: Callable[[str], T] = disabled_func("loader", "dumper")
    _dumper: Callable[[T, str], Any] = disabled_func("dumper", "loader")

    def dump(self, obj: T):
        """Write the supplied object using the configured dumper."""
        self._dumper(obj, self)

    def load(self):
        """Read the file using the configured loader and return its result."""
        return self._loader(self)  # type: ignore[arg-type]


class TypedFile(_TypedFile, Generic[T]):
    """
    A file path associated with a specific type for structured serialization.

    Behaves as a string file path, with ``dump()`` and ``load()`` methods
    that construct instances of the associated type when writing and reading.
    Custom callable dumpers and loaders replace these methods entirely.

    Formats that serialize mappings require instances to be convertible to
    dictionaries (see :func:`typed_to_dict`).
    Type annotations serve as documentation only; no runtime type checking is performed.
    """

    type_: Type[T]

    def __new__(cls, value: str, type_: Type[T]):
        self = super().__new__(cls, value)
        self.type_ = type_
        return self

    def dump(self, *args, **kwargs):
        """Construct ``type_(*args, **kwargs)`` and write it using the dumper."""
        obj = self.type_(*args, **kwargs)
        self._dumper(obj, self)

    def load(self):
        """Read a mapping using the loader and return ``type_(**data)``."""
        return self.type_(**self._loader(self))  # type: ignore[arg-type]


@overload
def typed_factory(
    type_: "None | str" = None,
    *,
    loader: "None | str | Callable[[str], T]" = None,
    dumper: "None | str | Callable[[T, str], Any]" = None,
) -> Callable[[str], _TypedFile[T]]: ...
@overload
def typed_factory(
    type_: "Type[T]",
    *,
    loader: "None | str | Callable[[str], T]" = None,
    dumper: "None | str | Callable[[T, str], Any]" = None,
) -> Callable[[str], TypedFile[T]]: ...
def typed_factory(  # type: ignore[reportInconsistentOverloads]
    type_: "None | str | Type[T]" = None,
    *,
    loader: "None | str | Callable[[str], T]" = None,
    dumper: "None | str | Callable[[T, str], Any]" = None,
):
    """
    Return a constructor that accepts a path and creates a string subclass.

    The subclass behaves as a plain string (file path) in all contexts,
    but provides ``dump()`` and ``load()`` methods to write and read
    instances of the associated type to and from disk.
    Constructing it does not perform I/O.

    ``type_`` controls object construction and default format selection:
    ---
    * ``Type`` (a class).
      ``load()`` constructs ``T(**data)`` from the decoded mapping and
      ``dump(*args, **kwargs)`` constructs ``T(*args, **kwargs)`` before writing.
    * ``str`` (A file format name), such as "json", "csv", or "json.gz",
      supplies defaults for any unspecified ``loader`` and ``dumper``
      via ``resolve_file_format`` (see below).
      In this case, supplying both ``loader`` and ``dumper`` is invalid.

    ``loader`` selects how to read the file:
    ---
    * ``str`` (A file format name)
      The decoded ``value`` is returned as ``T(**value)``
      when ``type_`` is a class, or directly otherwise.
    * ``Callable[[str], T]`` (receives the path),
      e.g., ``pd.read_csv`` or a class method ``T.load``.
      It overrides ``load()`` to return its result directly.

    ``dumper`` selects how to write the file:
    ---
    * ``str`` (A file format name).
      Formats such as JSON, YAML, and TOML convert the object to a dict.
      convert ``obj: T`` using ``obj.asdict()``
      if ``T`` is not dataclass, NamedTuple, or a Pydantic model.
    * ``Callable[[T, str], Any]`` (save ``T`` to path),
      e.g., ``pd.DataFrame.to_csv``.
      It overrides ``dump()`` to accept an already constructed object.

    If any of ``loader`` or ``dumper`` is not set:
    ---
    * ``(type_: str)`` (only ``type_`` supplied a file format name)
      equivalent to ``(loader=type_, dumper=type_)``.
    * ``(type_: str, loader: str | Callable)``
      and ``(type_: str, dumper: str | Callable)``
      equivalent to ``(loader=loader, dumper=type_)``
      and ``(loader=type_, dumper=dumper)``, respectively.
    * ``(type_: Type[T])`` (only ``type_`` supplied a class)
      both are inferred from the path suffix.
    * ``(type_: Type[T] | None, dumper: str)``
      (``type_`` is a class or is omitted, with a file format name for ``dumper``)
      equivalent to ``(type_=type_, loader=dumper, dumper=dumper)``.
    * ``(type_: Type[T] | None, dumper: Callable)``
      (``type_`` is a class or is omitted, with a callable for ``dumper``)
      dumper is replaced entirely, so loader cannot be inferred automatically.
      a ``ValueError`` is raised if ``load()`` is called
    * ``(type_: Type[T] | None, loader: Callable | str)``
      loader is replaced entirely, so dumper cannot be inferred automatically.
      a ``ValueError`` is raised if ``dump(*args, **kwargs)`` is called
    """
    typed_ = _TypedFile
    if type_ is None:
        if not (loader or dumper):
            raise ValueError(
                "At least one of type_, loader, or dumper must be supplied"
            )
    elif isinstance(type_, str):
        if not loader:
            loader = type_
        elif dumper:
            raise ValueError(
                f"Declared file format '{type_}' is completely overridden"
                " by supplied loader and dumper,"
                " please omit it or supply only one of loader or dumper"
            )
        dumper = dumper or type_
    else:
        typed_ = lambda file: TypedFile(file, type_)  # type: ignore[assignment]

    def _(file: str):
        typed_file: _TypedFile[T] = typed_(file)
        if dumper:
            if isinstance(dumper, str):
                default_loader, typed_file._dumper = resolve_file_format(dumper)
                if not loader:
                    typed_file._loader = default_loader
            else:
                typed_file.dump = lambda obj: dumper(obj, typed_file)  # type: ignore[method-assign, misc]

        if loader:
            if isinstance(loader, str):
                if loader == dumper:
                    typed_file._loader = default_loader
                else:
                    typed_file._loader = resolve_file_format(loader)[0]
            else:
                typed_file.load = lambda: loader(typed_file)  # type: ignore[method-assign]
        elif not dumper:
            typed_file._loader, typed_file._dumper = resolve_file_format(file)

        return typed_file

    return _


if sys.version_info < (3, 11):

    def typed_to_dict(obj):
        if isinstance(obj, tuple):
            if hasattr(obj, "_asdict"):
                return obj._asdict()
        elif hasattr(obj, "asdict"):
            return obj.asdict()
        raise NotImplementedError(f"Cannot convert {type(obj)} to dict")

    def resolve_file_format(
        file: str,
    ) -> Tuple[Callable[[str], Any], Callable[[Any, str], None]]:
        suffix = file.rsplit(".", 1)[-1]
        open_: Callable[[str, str], Any]
        if suffix == "gz":
            import gzip

            open_ = lambda f, mode: gzip.open(f, f"{mode}b")
        else:
            open_ = open
        if open_ is not open:
            suffix = file.rsplit(".", 2)[-2]
        if suffix in ("yaml", "yml"):
            import yaml

            return (
                lambda f: yaml.safe_load(open_(f, "r")),
                lambda obj, f: yaml.safe_dump(typed_to_dict(obj), stream=open_(f, "w")),  # type: ignore
            )
        elif suffix == "csv":
            import pandas as pd

            return (
                lambda f: pd.read_csv(f),
                lambda obj, f: obj.to_csv(f, index=False),  # type: ignore
            )
        elif suffix == "tsv":
            import pandas as pd

            return (
                lambda f: pd.read_csv(f, sep="\t"),
                lambda obj, f: obj.to_csv(f, sep="\t", index=False),  # type: ignore
            )
        raise NotImplementedError(f"Unsupported format[{suffix}] of file[{file}]")

else:
    from snakemake.io.resolve_file_format import resolve_file_format
