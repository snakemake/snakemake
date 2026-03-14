import os
from pathlib import Path
from typing import Callable, Optional, Union

from snakemake.common import get_function_params, overwrite_function_params
from snakemake.io import is_callable


def subpath(
    path_or_func: Union[Callable, str, Path],
    strip_suffix: Optional[str] = None,
    basename: bool = False,
    parent: bool = False,
    ancestor: Optional[int] = None,
):
    """Return the subpath of a given path.

    Args:
        path_or_func: A string, pathlib.Path, or a function returning a string or pathlib.Path.
        strip_suffix: If given, strip the suffix from the path.
        basename: If True, return the basename of the path (cannot be used together with parent or ancestor).
        parent: If True, return the parent directory of the path (cannot be used together with ancestor).
        ancestor: If given, return the ancestor directory of the path (cannot be used together with parent).
    """

    # TODO typecheck arguments, maybe find a decorator for that purpose
    def do(path):
        if isinstance(path, Path):
            path = str(path)
        if not isinstance(path, str):
            raise ValueError(
                "Value passed to subpath "
                "must be a single string or pathlib.Path (or a function returning those). "
                f"Obtained value: {repr(path)}"
            )
        if strip_suffix is not None:
            if not path.endswith(strip_suffix):
                raise ValueError(
                    f"Path {path} does not end with the specified suffix {strip_suffix}"
                )
            path = path[: -len(strip_suffix)]
        if basename:
            if parent or ancestor is not None:
                raise ValueError(
                    "basename cannot be used together with parent or ancestor"
                )
            path = os.path.basename(path)
        elif parent:
            if ancestor is not None:
                raise ValueError("parent cannot be used together with ancestor")
            path = os.path.dirname(path)
            if path == "":
                path = "."
        elif ancestor is not None:
            if ancestor < 1:
                raise ValueError("ancestor must be greater than 0")
            for _ in range(ancestor):
                path = os.path.dirname(path)
                if path == "":
                    path = "."
                    break
        return path

    if is_callable(path_or_func):
        params = get_function_params(path_or_func)

        def inner(wildcards, **args):
            value = path_or_func(wildcards, **args)
            return do(value)

        overwrite_function_params(inner, params)
        return inner
    else:
        return do(path_or_func)
