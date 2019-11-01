from snakemake.exceptions import IncompleteCheckpointException, WorkflowError
from snakemake.io import checkpoint_target


class Checkpoints:
    """ A namespace for checkpoints so that they can be accessed via dot notation. """

    def __init__(self):
        self.future_output = None

    def register(self, rule):
        setattr(self, rule.name, Checkpoint(rule, self))


class Checkpoint:
    __slots__ = ["rule", "checkpoints"]

    def __init__(self, rule, checkpoints):
        self.rule = rule
        self.checkpoints = checkpoints

    def get(self, **wildcards):
        missing = self.rule.wildcard_names.difference(wildcards.keys())
        if missing:
            raise WorkflowError(
                "Missing wildcard values for {}".format(", ".join(missing))
            )

        output, _ = self.rule.expand_output(wildcards)
        if self.checkpoints.future_output is None or any(
            (not f.exists or f in self.checkpoints.future_output) for f in output
        ):
            raise IncompleteCheckpointException(self.rule, checkpoint_target(output[0]))

        return CheckpointJob(self.rule, output)


class CheckpointJob:
    __slots__ = ["rule", "output"]

    def __init__(self, rule, output):
        self.output = output
        self.rule = rule
