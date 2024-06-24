from abc import ABC, abstractmethod
from collections import namedtuple
from collections.abc import Mapping, Callable
from functools import partial
import inspect
import os
import re
from typing import List, Optional, Union

from snakemake.common import async_run
import snakemake.io
import snakemake.utils
from snakemake.exceptions import LookupError
from snakemake_interface_common.exceptions import WorkflowError


def evaluate(expr: str):
    """Evaluate a python expression while replacing any wildcards given as
    {wildcardname} with the wildcard value represented as a string."""

    def inner(wildcards):
        return eval(
            expr.format(**{w: repr(v) for w, v in wildcards.items()}), globals()
        )

    return inner
