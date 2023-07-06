from collections.abc import Mapping, Callable
from typing import Optional, Union

import snakemake.io


def lookup(query: Optional[str] = None, dpath: Optional[str] = None, within=None):
    if within is None:
        raise ValueError(
            "Must provide a dataframe, series, or mapping to search within."
        )
    if query is not None:
        if isinstance(within, Mapping):
            raise ValueError(
                "Query parameter can only be used with pandas DataFrame or Series objects."
            )

        import pandas as pd

        if isinstance(within, pd.Series):
            within = within.to_frame()

        return within.query(query)
    elif dpath is not None:
        if not isinstance(within, Mapping):
            raise ValueError(
                "Dpath parameter can only be used with pandas DataFrame or Series objects."
            )
        import dpath

        return dpath.values(within, dpath)
    else:
        raise ValueError("Must provide either a query or dpath parameter.")


def evaluate(query: str):
    def inner(wildcards):
        return eval(query, wildcards)
    return inner


def branch(
    condition: Callable,
    then: Optional[Union[str, list[str], Callable]] = None,
    otherwise: Optional[Union[str, list[str], Callable]] = None,
):

    def inner(wildcards):
        def handle_callable(value):
            if isinstance(value, Callable):
                return value(wildcards)
            else:
                return value

        if condition(wildcards):
            return handle_callable(then)
        else:
            return handle_callable(otherwise)

    return inner


collect = snakemake.io.expand