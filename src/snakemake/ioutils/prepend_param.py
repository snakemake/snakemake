#!/usr/bin/env python3

from pathlib import Path
from typing import Callable, List, Optional, Union

from snakemake.common import get_function_params, overwrite_function_params
from snakemake.io import is_callable
from snakemake.ioutils.rule_items_proxy import RuleItemProxy


def prepend_param(
    prefix: str,
    paths_or_func: Union[Callable, str, Path, List[Union[str, Path]]],
    space: bool = True,
):
    """
    Prepend each filename in `paths_or_func` with `prefix`,
    optionally with a space (if `space` is `True`),
    and return a single combined string.
    This allows easier use of
    tools requiring multiple filenames to be prepended with a keyword.
    """

    def do_single(path):
        if isinstance(path, Path):
            path = str(path)
        if not isinstance(path, str):
            raise ValueError(
                "Values passed to prepend "
                "must be a single string or pathlib.Path (or a function returning those). "
                f"Obtained value: {repr(path)}"
            )
        return f"{prefix}{' ' if space else ''}{path}"

    def do(paths):
        if isinstance(paths, RuleItemProxy):

            def inner(wildcards, **kwargs):
                return do(kwargs[paths.name])

            return inner

        if not isinstance(paths, list):
            paths = [paths]

        return " ".join(map(do_single, paths))

    if is_callable(paths_or_func):
        params = get_function_params(paths_or_func)

        def inner(wildcards, **args):
            value = paths_or_func(wildcards, **args)
            return do(value)

        overwrite_function_params(inner, params)
        return inner
    else:
        return do(paths_or_func)
