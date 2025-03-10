from snakemake.logging import format_dict
from snakemake_interface_common.exceptions import WorkflowError


def evaluate(expr: str):
    """Evaluate a python expression while replacing any wildcards given as
    {wildcardname} with the wildcard value represented as a string."""

    def inner(wildcards):
        formatted = expr.format(**{w: repr(v) for w, v in wildcards.items()})
        try:
            return eval(formatted, globals())
        except Exception as e:
            raise WorkflowError(
                f"Failed to evaluate expression {expr} with wildcards {format_dict(wildcards)}. "
                f"Formatted expression: {formatted}",
                e,
            )

    return inner
