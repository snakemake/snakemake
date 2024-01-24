from collections import namedtuple
from collections.abc import Mapping, Callable
from typing import List, Optional, Union

import snakemake.io
import snakemake.utils
from snakemake_interface_common.exceptions import WorkflowError


def lookup(
    dpath: Optional[str] = None,
    query: Optional[str] = None,
    cols: Optional[List[str]] = None,
    is_nrows: Optional[int] = None,
    within=None,
):
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
    Further, it is possible to constrain the output to a list of columns, e.g.
    `lookup(query="index.loc[{sample}]", within=samples, cols=["somecolumn"])`.

    In case of a pandas series, the series is converted into a dataframe via
    Series.to_frame() and the same logic as for a dataframe is applied.

    In case of a python mapping, the dpath parameter is passed to dpath.values()
    (see https://github.com/dpath-maintainers/dpath-python).

    Query, dpath and cols may contain wildcards (e.g. {sample}).
    In that case, this function returns a Snakemake input function which takes
    wildcards as its only argument and will be evaluated by Snakemake
    once the wildcard values are known.
    """
    if within is None:
        raise ValueError(
            "Must provide a dataframe, series, or mapping to search within."
        )

    def handle_wildcards(expression, func, cols=None, is_nrows=None):
        if snakemake.io.contains_wildcard(expression) or (
            cols is not None
            and any(snakemake.io.contains_wildcard(col) for col in cols)
        ):

            def inner(wildcards):
                resolved_expression = snakemake.utils.format(expression, **wildcards)
                if cols is not None:
                    resolved_cols = [snakemake.utils.format(col, **wildcards) for col in cols]
                    return func(resolved_expression, cols=resolved_cols, is_nrows=is_nrows)
                else:
                    return func(resolved_expression)

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

        def do_query(query, cols=None, is_nrows=None):
            try:
                res = within.query(query)
            except Exception as e:
                raise WorkflowError(f"Error in lookup function", e)
            if isinstance(res, pd.Series):
                if is_nrows is not None:
                    return is_nrows == 1
                if cols is not None:
                    res = res[cols]
                # convert series into named tuple with index as attribute
                thetuple = namedtuple("Row", ["index"] + res.index.tolist())(**res)
                return thetuple
            else:
                if is_nrows is not None:
                    return is_nrows == len(res)
                if cols is not None:
                    res = res[cols]
                res = list(res.itertuples())
                if len(res) == 1:
                    # just return the item if it is only one
                    return res[0]
                return res

        return handle_wildcards(query, do_query, cols=cols, is_nrows=is_nrows)

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
