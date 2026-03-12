#!/usr/bin/env python3

import re


def format_python_module(filename: str) -> str:
    """
    Given the path to a `.py` file,
    return the name of the corresponding module that can be passed to `python -m`.
    For example,

    >>> python_module("package/subpackage/module.py")
    'package.subpackage.module'
    """

    if not filename.endswith(".py"):
        raise ValueError("Only .py files may be run as Python modules.")

    components = filename[:-3].split("/")
    if any(
        re.match("[A-Za-z_][A-Za-z0-9_]*", component) is None
        for component in components
    ):
        raise ValueError(f"{filename} does not translate to a valid Python name")

    return ".".join(components)


def python_module(key: str = "script"):
    """
    Returns a helper function to translate an input `.py` filename into a module name.
    By default,
    this uses `input.script`;
    the `key` argument allows another input file to be used.
    """
    return lambda wildcards, input: format_python_module(input[key])
