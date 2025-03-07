from typing import TYPE_CHECKING

from snakemake.exceptions import IncompleteCheckpointException, WorkflowError
from snakemake.io import checkpoint_target

if TYPE_CHECKING:
    from snakemake.rules import Rule


class Checkpoints:
    """A namespace for checkpoints so that they can be accessed via dot notation."""

    def __init__(self):
        self.created_output = None

    def register(self, rule, fallback_name=None):
        checkpoint = Checkpoint(rule, self)
        if fallback_name:
            setattr(self, fallback_name, checkpoint)
        setattr(self, rule.name, Checkpoint(rule, self))


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

        if self.checkpoints.created_output is not None:
            if set(output) <= set(self.checkpoints.created_output):
                return CheckpointJob(self.rule, output)

        raise IncompleteCheckpointException(self.rule, checkpoint_target(output[0]))


class CheckpointJob:
    __slots__ = ["rule", "output"]

    def __init__(self, rule: "Rule", output):
        self.output = output
        self.rule = rule
