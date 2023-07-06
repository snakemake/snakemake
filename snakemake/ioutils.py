from collections.abc import Mapping, Callable
from typing import Optional, Union

import snakemake.io
import snakemake.utils


def lookup(query: Optional[str] = None, dpath: Optional[str] = None, within=None):
    """Lookup values in a pandas dataframe, series, or python mapping (e.g. dict).

    In case of a pandas dataframe (see https://pandas.pydata.org),
    the query parameter is passed to DataFrame.query()
    and the result is returned as an iterator over the rows of the dataframe.
    Each row is represented as a named tuple with the column names as attributes.

    In case of a pandas series, the series is converted into a dataframe via
    Series.to_frame() and the same logic as for a dataframe is applied.

    In case of a python mapping, the dpath parameter is passed to dpath.values()
    (see https://github.com/dpath-maintainers/dpath-python).

    Both query and dpath may contain wildcards (e.g. {sample}).
    In that case, this function returns a Snakemake input function which takes
    wildcards as its only argument and will be evaluated by Snakemake
    once the wildcard values are known.
    """
    if within is None:
        raise ValueError(
            "Must provide a dataframe, series, or mapping to search within."
        )

    def handle_wildcards(expression, func):
        if snakemake.io.contains_wildcard(expression):

            def inner(wildcards):
                return func(snakemake.utils.format(expression, **wildcards))

            return inner
        else:
            return func(expression)

    if query is not None:
        if isinstance(within, Mapping):
            raise ValueError(
                "Query parameter can only be used with pandas DataFrame or Series objects."
            )

        import pandas as pd

        if isinstance(within, pd.Series):
            within = within.to_frame()

        def do_query(query):
            return within.query(query).itertuples()

        return handle_wildcards(query, do_query)

    elif dpath is not None:
        if not isinstance(within, Mapping):
            raise ValueError(
                "Dpath parameter can only be used with pandas DataFrame or Series objects."
            )
        import dpath as dp

        def do_dpath(dpath):
            return dp.values(within, dpath)

        return handle_wildcards(dpath, do_dpath)
    else:
        raise ValueError("Must provide either a query or dpath parameter.")


def evaluate(query: str):
    """Evaluate a python expression while replacing any wildcards given as
    {wildcardname} with the wildcard value represented as a string."""

    def inner(wildcards):
        return eval(query, {w: repr(v) for w, v in wildcards.items()})

    return inner


def branch(
    condition: Callable,
    then: Optional[Union[str, list[str], Callable]] = None,
    otherwise: Optional[Union[str, list[str], Callable]] = None,
):
    """Branch based on a condition that is provided as a function pointer (i.e. a Callable).
    If the condition is true, return the value
    of the then parameter. Otherwise, return the value of the otherwise parameter.

    The given condition function has to take wildcards as its only parameter.
    This function returns a derived input function that will be evaluated once the
    wildcards are known.
    """

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


# Alias for expand that provides a more intuitive name for the use case of
# collecting files from previous jobs.
collect = snakemake.io.expand


def register_in_globals(_globals):
    _globals.update(
        {
            "lookup": lookup,
            "evaluate": evaluate,
            "branch": branch,
            "collect": collect,
        }
    )
