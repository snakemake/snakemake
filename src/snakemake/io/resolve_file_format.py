import dataclasses
from typing import Callable, Any, Tuple


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


def _compressed_mode(mode: str) -> str:
    return mode if "b" in mode else f"{mode}t"


def resolve_file_format(
    file: str,
) -> Tuple[Callable[[str], Any], Callable[[Any, str], None]]:
    suffix = file.rsplit(".", 1)[-1]
    open_: Callable[[str, str], Any]
    match suffix:
        case "gz":
            import gzip

            open_ = lambda f, mode: gzip.open(f, _compressed_mode(mode))
        case "bz2":
            import bz2

            open_ = lambda f, mode: bz2.open(f, _compressed_mode(mode))
        case "xz":
            import lzma

            open_ = lambda f, mode: lzma.open(f, _compressed_mode(mode))
        case _:
            open_ = open
    if open_ is not open:
        suffix = file.rsplit(".", 2)[-2]
    match suffix:
        case "json":
            import json

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
                    _,
                )
        case "pkl" | "pickle":
            import pickle

            return (
                lambda f: pickle.load(open_(f, "rb")),  # type: ignore[arg-type]
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
