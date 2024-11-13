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


class WildcardHandlerBase(ABC):
    fmt_regex = re.compile(r"\{(?P<stmt>[^\{][^\{\}]+)\}[^\}]")

    def __init__(self, func, **namespace):
        self.func = func
        self.namespace = namespace

    def needs_wildcards(self, expression):
        return callable(expression) or any(
            name not in self.namespace
            for name in snakemake.io.get_wildcard_names(expression)
        )

    @abstractmethod
    def apply_func(self, expression, namespace=None): ...

    def handle(self, expression):
        if self.needs_wildcards(expression) or any(
            callable(value) for value in self.namespace.values()
        ):

            def inner(wildcards):
                if self.namespace:
                    # add wildcard values to namespace
                    # do not override namespace
                    # (as it has been chosen explicitly by the dev)
                    namespace = dict(self.namespace)
                    for name, value in list(namespace.items()):
                        # resolve callables given in namespace
                        if callable(value):
                            namespace[name] = value(wildcards)
                    for name, value in wildcards.items():
                        if name not in namespace:
                            namespace[name] = value
                else:
                    namespace = wildcards
                if callable(expression):
                    resolved_expression = expression(wildcards)
                else:
                    resolved_expression = expression
                resolved_expression = snakemake.utils.format(
                    resolved_expression, **namespace
                )
                return self.apply_func(resolved_expression, namespace)

            return inner
        else:
            if self.namespace:
                resolved_expression = snakemake.utils.format(
                    expression, **self.namespace
                )
            else:
                resolved_expression = expression
            return self.apply_func(resolved_expression, self.namespace)


class DpathWildcardHandler(WildcardHandlerBase):
    def apply_func(self, expression, namespace=None):
        return self.func(expression)


class QueryWildcardHandler(WildcardHandlerBase):
    def __init__(self, func, cols=None, is_nrows=None, **namespace):
        super().__init__(func, **namespace)
        self.cols = cols
        self.is_nrows = is_nrows

    def needs_wildcards(self, expression):
        if super().needs_wildcards(expression):
            return True
        if self.cols is None:
            return False
        if isinstance(self.cols, list):
            return any(
                super(QueryWildcardHandler, self).needs_wildcards(col)
                for col in self.cols
            )
        else:
            return super().needs_wildcards(self.cols)

    def apply_func(self, expression, namespace=None):
        cols = self.cols
        if self.cols is not None and namespace is not None:
            if isinstance(self.cols, list):
                cols = [snakemake.utils.format(col, **namespace) for col in self.cols]
            else:
                cols = snakemake.utils.format(self.cols, **namespace)
        return self.func(expression, cols=cols, is_nrows=self.is_nrows)


NODEFAULT = object()


def lookup(
    dpath: Optional[str] = None,
    query: Optional[str] = None,
    cols: Optional[Union[List[str], str]] = None,
    is_nrows: Optional[int] = None,
    within=None,
    default=NODEFAULT,
    **namespace,
):
    """Lookup values in a pandas dataframe, series, or python mapping (e.g. dict).

    Required argument ``within`` should be a pandas dataframe or series (in which
    case use ``query``, and optionally ``cols`` and ``is_nrows``), or a Python
    mapping like a dict (in which case use the ``dpath`` argument is used).

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
    ``branch(lookup(query="sample == '{sample}' & someval > 2", within=samples), then="foo", otherwise="bar")``.
    In case your dataframe has an index, you can also access the index within the
    query, e.g. for faster, constant time lookups: ``lookup(query="index.loc[{sample}]", within=samples)``.
    Further, it is possible to constrain the output to a list of columns, e.g.
    ``lookup(query="sample == '{sample}'", within=samples, cols=["somecolumn"])`` or to
    a single column, e.g.
    ``lookup(query="sample == '{sample}'", within=samples, cols="somecolumn")``.
    In the latter case, just a list of items in that column is returned.
    Finally, if the integer argument ``is_nrows`` is used, this returns true
    if there are that many rows in the query results, false otherwise.

    In case of a pandas series, the series is converted into a dataframe via
    Series.to_frame() and the same logic as for a dataframe is applied.

    In case of a python mapping, the ``dpath`` parameter is passed to
    ``dpath.values()`` (see https://github.com/dpath-maintainers/dpath-python),
    and the ``query``, ``cols``, and ``is_nrows`` arguments are ignored. If the
    dpath is not found, a ``LookupError`` is raised, unless a default fallback
    value is provided via the ``default`` argument.

    Query, dpath and cols may contain wildcards (e.g. {sample}).
    In that case, this function returns a Snakemake input function which takes
    wildcards as its only argument and will be evaluated by Snakemake
    once the wildcard values are known.

    In addition to wildcard values, dpath, query and cols may refer via the same syntax
    to auxiliary namespace arguments given to the lookup function, e.g.
    ``lookup(query="cell_type == '{sample.cell_type}'", within=samples, sample=lookup("sample == '{sample}'", within=samples))``
    This way, one can e.g. pass additional variables or chain lookups into more complex queries.
    """
    error = partial(LookupError, query=query, dpath=dpath)

    if within is None:
        raise error(
            msg="Must provide a dataframe, series, or mapping to search within."
        )
    if cols is not None and not isinstance(cols, (str, list)):
        raise error(msg="The cols argument has to be either a str or a list of str.")
    if is_nrows is not None and not isinstance(is_nrows, int):
        raise error(msg="The is_nrows argument has to be an int.")

    if query is not None:
        if isinstance(within, Mapping):
            raise error(
                msg=f"Query parameter can only be used with pandas DataFrame or Series objects."
            )

        import pandas as pd

        if isinstance(within, pd.Series):
            within = within.to_frame()

        def do_query(query, cols=None, is_nrows=None):
            try:
                res = within.query(query)
            except Exception as e:
                raise LookupError(query=query, exc=e)

            if is_nrows is not None:
                return is_nrows == len(res)
            if cols is not None:
                res = res[cols]
                if not isinstance(cols, list):
                    # single column select, just return a list of values
                    return res.to_list()
            res = list(res.itertuples(index=cols is None))
            if len(res) == 1:
                # just return the item if it is only one
                return res[0]
            return res

        return QueryWildcardHandler(
            do_query, cols=cols, is_nrows=is_nrows, **namespace
        ).handle(query)

    elif dpath is not None:
        if not isinstance(within, Mapping):
            raise error(
                msg="Dpath parameter can only be used with python mapping (e.g. dict)."
            )
        import dpath as dp

        def do_dpath(dpath):
            try:
                return dp.get(within, dpath)
            except ValueError:
                return dp.values(within, dpath)
            except KeyError:
                if default is not NODEFAULT:
                    return default
                raise LookupError(dpath=dpath, msg="Dpath not found.")

        return DpathWildcardHandler(do_dpath, **namespace).handle(dpath)
    else:
        raise error("Must provide either a query or dpath parameter.")


def evaluate(expr: str):
    """Evaluate a python expression while replacing any wildcards given as
    {wildcardname} with the wildcard value represented as a string."""

    def inner(wildcards):
        formatted = expr.format(**{w: repr(v) for w, v in wildcards.items()})
        try:
            return eval(formatted, globals())
        except Exception as e:
            raise WorkflowError(
                f"Failed to evaluate expression {expr} with wildcards {wildcards}. Formatted expression: {formatted}",
                e,
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


def exists(path):
    """Return True if the given file or directory exists.

    This function considers any storage arguments given to Snakemake.
    """
    func_context = inspect.currentframe().f_back.f_locals
    func_context_global = inspect.currentframe().f_back.f_globals

    workflow = func_context.get("workflow") or func_context_global.get("workflow")

    if workflow is None:
        raise WorkflowError(
            "The exists function can only be used within a Snakemake workflow "
            "(the global variable 'workflow' has to be present)."
        )

    path = workflow.modifier.path_modifier.apply_default_storage(path)
    if snakemake.io.is_flagged(path, "storage_object"):
        return async_run(path.flags["storage_object"].managed_exists())
    else:
        return os.path.exists(path)


def register_in_globals(_globals):
    _globals.update(
        {
            "lookup": lookup,
            "evaluate": evaluate,
            "branch": branch,
            "collect": collect,
            "exists": exists,
        }
    )
