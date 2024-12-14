from collections.abc import Mapping, Callable
from typing import Optional, Union


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
        try:
            selected_case = cases[res]
        except KeyError:
            raise KeyError(f"Key {res} not found in given cases of branch function")
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
