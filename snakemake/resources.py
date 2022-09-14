import re


class DefaultResources:
    def __init__(self, args=None):
        self.args = args

        def fallback(val):
            def callable(wildcards, input, attempt, threads, rulename):
                try:
                    value = eval(
                        val, {"input": input, "attempt": attempt, "threads": threads}
                    )
                # Triggers for string arguments like n1-standard-4
                except NameError:
                    return val
                return value

            return callable

        self.parsed = dict(_cores=1, _nodes=1)
        if self.args is not None:
            self.parsed.update(parse_resources(args, fallback=fallback))


def parse_resources(resources_args, fallback=None):
    """Parse resources from args."""
    resources = dict()
    if resources_args is not None:
        valid = re.compile(r"[a-zA-Z_]\w*$")
        for res in resources_args:
            try:
                res, val = res.split("=")
            except ValueError:
                raise ValueError("Resources have to be defined as name=value pairs.")
            if not valid.match(res):
                raise ValueError(
                    "Resource definition must start with a valid identifier."
                )
            try:
                val = int(val)
            except ValueError:
                if fallback is not None:
                    val = fallback(val)
                else:
                    raise ValueError(
                        "Resource definiton must contain an integer after the identifier."
                    )
            if res == "_cores":
                raise ValueError(
                    "Resource _cores is already defined internally. Use a different name."
                )
            resources[res] = val
    return resources
