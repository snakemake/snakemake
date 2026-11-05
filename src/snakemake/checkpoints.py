from typing import TYPE_CHECKING
from snakemake.exceptions import IncompleteCheckpointException, WorkflowError
from snakemake.io import checkpoint_target
from snakemake.logging import logger

if TYPE_CHECKING:
    from snakemake.rules import Rule


class Checkpoints:
    """A singleton object in a workflow.

    Created_output can be accessed by checkpoint rules in or out modules.
    This never go into snakefile, so no rules name will be set to it.
    """

    def __init__(self):
        self._created_output = set()
        self._rules: dict[str, "Rule"] = {}

    @property
    def created_output(self):
        return self._created_output

    def register(self, rule: "Rule", name: "str | None" = None):
        self._rules[name or rule.name] = rule

    def _resolve(self, module_name: "str | None", name: str):
        try:
            rule = self._rules[name]
        except KeyError:
            raise AttributeError(name)
        return Checkpoint(rule, self)

    def _spawn_new_namespace(
        self, module_name: "str | None" = None, parent: "CheckpointsProxy | None" = None
    ):
        """Make a new namespace for checkpoints in the module."""
        return CheckpointsProxy(parent or self, module_name)


class CheckpointsProxy:
    """A view on the workflow checkpoints for one module namespace.

    The proxy stores no checkpoints itself. It only maps local checkpoint
    names to names in its parent namespace. Resolution walks up the parent
    chain until the workflow's true checkpoint registry is reached.
    """

    def __init__(
        self, parent: "Checkpoints | CheckpointsProxy", module_name: "str | None" = None
    ):
        self._parent = parent
        self._module_name = module_name
        self._child_refs: "dict[tuple[str | None, str], str]" = {}

    def _register(self, module_name: "str | None", orig_name: str, name: str):
        key = module_name, orig_name
        if key in self._child_refs:
            if self._child_refs[key] != name:
                logger.warning(
                    f"Checkpoint {orig_name!r} of module {module_name!r} is already "
                    f"mapped to {self._child_refs[key]!r} in the parent namespace, "
                    f"this mapping will not overwrite the existing mapping."
                )
            return
        self._child_refs[key] = name

    def _resolve(self, module_name: "str | None", name: str) -> "Checkpoint":
        key = module_name, name
        try:
            local_name = self._child_refs[key]
        except KeyError:
            raise AttributeError(
                f"Checkpoint {name!r} of module {module_name!r} is not mapped "
                "in the parent checkpoint namespace."
            )
        try:
            return self._parent._resolve(self._module_name, local_name)
        except AttributeError as err:
            raise WorkflowError(
                f"Checkpoint {name!r} points to missing checkpoint "
                f"{local_name!r} in the parent namespace."
            ) from err

    def get(self, name: str):
        return self._parent._resolve(self._module_name, name)

    def as_namespace(self):
        return type(
            "Checkpoints", (), {"__getattr__": lambda _, name: self.get(name)}
        )()


class Checkpoint:
    __slots__ = ["rule", "checkpoints"]

    def __init__(self, rule: "Rule", checkpoints: Checkpoints):
        self.rule = rule
        self.checkpoints = checkpoints

    def get(self, **wildcards):
        missing = self.rule.wildcard_names.difference(wildcards.keys())
        if missing:
            raise WorkflowError(
                "Missing wildcard values for {}".format(", ".join(missing))
            )

        output, _ = self.rule.expand_output(wildcards)
        if self.checkpoints.created_output:
            missing_output = set(output) - set(self.checkpoints.created_output)
            if not missing_output:
                return CheckpointJob(self.rule, output)
            else:
                logger.debug(
                    f"Missing checkpoint output for {self.rule.name} "
                    f"(wildcards: {wildcards}): {','.join(missing_output)} of {','.join(output)}"
                )

        raise IncompleteCheckpointException(self.rule, checkpoint_target(output[0]))


class CheckpointJob:
    __slots__ = ["rule", "output"]

    def __init__(self, rule: "Rule", output):
        self.output = output
        self.rule = rule
