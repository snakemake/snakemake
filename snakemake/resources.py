import re
import tempfile


class DefaultResources:
    defaults = {
        "mem_mb": "max(2*input.size_mb, 1000)",
        "disk_mb": "max(2*input.size_mb, 1000)",
        "tmpdir": "system_tmpdir",
    }

    @classmethod
    def decode_arg(cls, arg):
        try:
            return arg.split("=")
        except ValueError:
            raise ValueError("Resources have to be defined as name=value pairs.")

    @classmethod
    def encode_arg(cls, name, value):
        return "{}={}".format(name, value)

    def __init__(self, args=None):
        if args is None:
            args = []
        self._args = dict(DefaultResources.defaults)
        self._args.update({name: value for name, value in map(self.decode_arg, args)})

        def fallback(val):
            def callable(wildcards, input, attempt, threads, rulename):
                try:
                    value = eval(
                        val,
                        {
                            "input": input,
                            "attempt": attempt,
                            "threads": threads,
                            "system_tmpdir": tempfile.gettempdir(),
                        },
                    )
                # Triggers for string arguments like n1-standard-4
                except NameError:
                    return val
                return value

            return callable

        self.parsed = dict(_cores=1, _nodes=1)
        self.parsed.update(parse_resources(self._args, fallback=fallback))

    @property
    def args(self):
        return [self.encode_arg(name, value) for name, value in self._args.items()]


def parse_resources(resources_args, fallback=None):
    """Parse resources from args."""
    resources = dict()
    if resources_args is not None:
        valid = re.compile(r"[a-zA-Z_]\w*$")
        for res, val in resources_args.items():
            if not valid.match(res):
                raise ValueError(
                    "Resource definition must start with a valid identifier, but found {}.".format(
                        res
                    )
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
