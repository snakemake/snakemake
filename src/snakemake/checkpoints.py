import threading
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

    @property
    def created_output(self):
        return self._created_output

    def spawn_new_namespace(self):
        """Make a new namespace for checkpoints in the module."""
        return CheckpointsProxy(self)


class CheckpointsProxy:
    """A namespace for checkpoints so that they can be accessed via dot notation.

    It will be created once a module is created,
    and different module will have different checkpoint namespace,
    but share a single created_output set.
    """

    def __init__(self, parent: Checkpoints):
        self._parent = parent
        self._rules = {}
        self._local = threading.local()

    def __getattr__(self, name):
        if name in self._rules:
            checkpoint = self._rules[name]
            if getattr(self._local, "cache", None):
                return CheckpointLater(checkpoint, self._local.cache[-1])
            return checkpoint
        raise WorkflowError(
            f"Checkpoint {name} is not defined in this workflow. "
            f"Available checkpoints: {', '.join(self._rules)}"
        )

    def register(self, rule: "Rule", fallback_name=None):
        self._rules[fallback_name or rule.name] = Checkpoint(rule, self._parent)

    def __enter__(self):
        if getattr(self._local, "cache", None) is None:
            self._local.cache = []
        # do not lost message with nested checkpoint
        self._local.cache.append([])
        return self

    def __exit__(self, exc_type, exc, tb):
        cache: "list[IncompleteCheckpointException]" = self._local.cache.pop()
        if exc_type is None and cache:
            e, *_e = cache
            for c in self._local.cache:
                _e.extend(c)
            e.extra = _e
            raise e
        return False


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
        return self.expect(output)

    def expect(self, output) -> "CheckpointJob":
        raise IncompleteCheckpointException(self.rule, checkpoint_target(output[0]))


class CheckpointLater(Checkpoint):
    def __init__(self, checkpoint: Checkpoint, exceptions: list):
        super().__init__(checkpoint.rule, checkpoint.checkpoints)
        self.exceptions = exceptions

    def expect(self, output):
        self.exceptions.append(
            IncompleteCheckpointException(self.rule, checkpoint_target(output[0]))
        )
        return CheckpointJob(self.rule, output)


class CheckpointJob:
    __slots__ = ["rule", "output"]

    def __init__(self, rule: "Rule", output):
        self.output = output
        self.rule = rule
