from collections import namedtuple
from collections.abc import Mapping, Callable
from typing import Optional, Union

import snakemake.io
import snakemake.utils
from snakemake_interface_common.exceptions import WorkflowError


def lookup(query: Optional[str] = None, dpath: Optional[str] = None, within=None):
    """Lookup values in a pandas dataframe, series, or python mapping (e.g. dict).

    In case of a pandas dataframe (see https://pandas.pydata.org),
    the query parameter is passed to DataFrame.query().
    If the query results in multiple rows, the result is returned as a list of
    named tuples with the column names as attributes.
    If the query results in a single row, the result is returned as a single
    named tuple with the column names as attributes.
    In both cases, the result can be used by the expand or collect function,
    e.g. `collect("results/{item.sample}.txt", sample=lookup(query="someval > 2", within=samples))`.
    Since the result, in any case, also evaluates to True if it is not empty
    when interpreted as a boolean by Python, it can also be used as a condition
    for the branch function, e.g.
    `branch(lookup(query="sample == {sample} & someval > 2", within=samples), then="foo", otherwise="bar")`.
    In case your dataframe has an index, you can also access the index within the
    query, e.g. for faster, constant time lookups: `lookup(query="index.loc[{sample}]", within=samples)`.

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
            try:
                res = within.query(query)
            except Exception as e:
                raise WorkflowError(f"Error in lookup function", e)
            if isinstance(res, pd.Series):
                # convert series into named tuple with index as attribute
                return namedtuple("Row", ["index"] + res.index.tolist())(**res)
            else:
                res = list(res.itertuples())
                if len(res) == 1:
                    # just return the item if it is only one
                    return res[0]
                return res

        return handle_wildcards(query, do_query)

    elif dpath is not None:
        if not isinstance(within, Mapping):
            raise ValueError(
                "Dpath parameter can only be used with pandas DataFrame or Series objects."
            )
        import dpath as dp

        def do_dpath(dpath):
            try:
                return dp.get(within, dpath)
            except ValueError:
                return dp.values(within, dpath)
            except KeyError as e:
                raise WorkflowError(
                    f"Error in lookup function: dpath {dpath} not found."
                )

        return handle_wildcards(dpath, do_dpath)
    else:
        raise ValueError("Must provide either a query or dpath parameter.")


def evaluate(expr: str):
    """Evaluate a python expression while replacing any wildcards given as
    {wildcardname} with the wildcard value represented as a string."""

    def inner(wildcards):
        return eval(
            expr.format(**{w: repr(v) for w, v in wildcards.items()}), globals()
        )

    return inner


def branch(
    condition: Union[Callable, bool],
    then: Optional[Union[str, list[str], Callable]] = None,
    otherwise: Optional[Union[str, list[str], Callable]] = None,
    cases: Optional[Mapping] = None,
):
    """Branch based on a condition that is provided as a function pointer (i.e. a Callable)
    or a value.

    If then and optionally otherwise are specified, do the following:
    If the condition is (or evaluates to) True, return the value
    of the then parameter. Otherwise, return the value of the otherwise parameter.

    If cases is specified, do the following:
    Retrieve the value of the cases mapping using the return value of the condition
    (if it is a function), or the condition value itself as a key.

    The given condition function has to take wildcards as its only parameter.
    Similarly, then, otherwise and the values of the cases mapping can be such functions.

    If any such function is given to any of those arguments, this function returns a derived
    input function that will be evaluated once the wildcards are known.
    """

    def convert_none(value):
        return value or []

    def handle_callable(value, wildcards):
        if isinstance(value, Callable):
            return convert_none(value(wildcards))
        else:
            return convert_none(value)

    def do_branch_then_otherwise(wildcards):
        if handle_callable(condition, wildcards):
            return handle_callable(then, wildcards)
        else:
            return handle_callable(otherwise, wildcards)

    def do_branch_cases(wildcards):
        res = handle_callable(condition, wildcards)
        selected_case = cases[res]
        return handle_callable(selected_case, wildcards)

    do_branch = do_branch_then_otherwise
    if cases is not None:
        if otherwise is not None or then is not None:
            raise ValueError("Cannot use cases together with then or otherwise.")
        do_branch = do_branch_cases

    if any(isinstance(value, Callable) for value in (condition, then, otherwise)):

        def inner(wildcards):
            return do_branch(wildcards)

        return inner
    else:
        return do_branch(None)


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
