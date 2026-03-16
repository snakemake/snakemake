#!/usr/bin/env python3

from pathlib import Path
import re
from typing import Callable, Union

from snakemake.common import get_function_params, overwrite_function_params
from snakemake.io import is_callable
from snakemake.ioutils.rule_items_proxy import rule_item_factory, RuleItemProxy


def format_python_module(filepath: Union[Path, str]) -> str:
    """
    Given the path to a `.py` file,
    return the name of the corresponding module that can be passed to `python -m`.
    For example,

    >>> format_python_module("package/subpackage/module.py")
    'package.subpackage.module'
    """

    filepath = Path(filepath)
    if not filepath.name.endswith(".py"):
        raise ValueError("Only .py files may be run as Python modules.")

    components = list(filepath.parts)

    # Strip `.py`
    components[-1] = filepath.name[:-3]

    if any(
        re.match("[A-Za-z_][A-Za-z0-9_]*$", component) is None
        for component in components
    ):
        raise ValueError(f"{filepath} does not translate to a valid Python name")

    return ".".join(components)


def as_py_module(
    path_or_func: Union[Callable, str, Path, RuleItemProxy] = rule_item_factory(
        "input"
    ).script,
):
    """
    Returns the Python module name associated with the `.py` file at a given path.
        path_or_func: A string, pathlib.Path, or a function returning a string or pathlib.Path.
                      Default: input.script
    """
    if is_callable(path_or_func):
        params = get_function_params(path_or_func)

        def inner(wildcards, **args):
            value = path_or_func(wildcards, **args)
            return as_py_module(value)

        overwrite_function_params(inner, params)
        return inner
    elif isinstance(path_or_func, RuleItemProxy):
        # This should only occur if `input` is the argument
        def inner(wildcards, **kwargs):
            # The rule item must resolve to a single path
            (filepath,) = kwargs[path_or_func.name]
            return format_python_module(filepath)

        return inner
    else:
        return format_python_module(path_or_func)

    return lambda wildcards, input: format_python_module(input[key])
