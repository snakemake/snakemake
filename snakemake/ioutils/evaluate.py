def evaluate(expr: str):
    """Evaluate a python expression while replacing any wildcards given as
    {wildcardname} with the wildcard value represented as a string."""

    def inner(wildcards):
        return eval(
            expr.format(**{w: repr(v) for w, v in wildcards.items()}), globals()
        )

    return inner
