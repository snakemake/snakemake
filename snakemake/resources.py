import re
import tempfile

from snakemake.exceptions import WorkflowError


class DefaultResources:
    defaults = {
        "mem_mb": "max(2*input.size_mb, 1000)",
        "disk_mb": "max(2*input.size_mb, 1000)",
        "tmpdir": "system_tmpdir",
    }

    bare_defaults = {
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

    def __init__(self, args=None, from_other=None, mode="full"):
        if mode == "full":
            self._args = dict(DefaultResources.defaults)
        elif mode == "bare":
            self._args = dict(DefaultResources.bare_defaults)
        else:
            raise ValueError("Unexpected mode for DefaultResources: {}".format(mode))

        if from_other is not None:
            self._args = dict(from_other._args)
            self.parsed = dict(from_other.parsed)
        else:
            if args is None:
                args = []

            self._args.update(
                {name: value for name, value in map(self.decode_arg, args)}
            )

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
                    except Exception as e:
                        if not (
                            isinstance(e, FileNotFoundError) and e.filename in input
                        ):
                            # Missing input files are handled by the caller
                            raise WorkflowError(
                                "Failed to evaluate default resources value "
                                "'{}'.\n"
                                "    String arguments may need additional "
                                "quoting. Ex: --default-resources "
                                "\"tmpdir='/home/user/tmp'\".".format(val),
                                e,
                            )
                        raise e
                    return value

                return callable

            self.parsed = dict(_cores=1, _nodes=1)
            self.parsed.update(parse_resources(self._args, fallback=fallback))

    def set_resource(self, name, value):
        self._args[name] = "{}".format(value)
        self.parsed[name] = value

    @property
    def args(self):
        return [self.encode_arg(name, value) for name, value in self._args.items()]

    def __bool__(self):
        return bool(self.parsed)


def parse_resources(resources_args, fallback=None):
    """Parse resources from args."""
    resources = dict()
    if resources_args is not None:
        valid = re.compile(r"[a-zA-Z_]\w*$")

        if isinstance(resources_args, list):
            resources_args = map(DefaultResources.decode_arg, resources_args)
        else:
            resources_args = resources_args.items()

        for res, val in resources_args:
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
