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


class CheckpointsProxy(Checkpoints):
    """A namespace for checkpoints so that they can be accessed via dot notation.

    It will be created once a module is created,
    and different module will have different checkpoint namespace,
    but share a single created_output set.
    """

    def __init__(self, parent: Checkpoints):
        self.parent = parent

    @property
    def created_output(self):
        return self.parent.created_output

    def register(self, rule: "Rule", fallback_name=None):
        checkpoint = Checkpoint(rule, self)
        if fallback_name:
            setattr(self, fallback_name, checkpoint)
        setattr(self, rule.name, checkpoint)


def dictproduct(listdict):
    """itertools.product() for dict of lists that yields a dict for each combination of values"""
    import itertools

    assert all(isinstance(v, list) for v in listdict.values())
    for x in itertools.product(*listdict.values()):
        yield dict(zip(listdict.keys(), x))


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

        listify = lambda x: (
            [
                x,
            ]
            if not isinstance(x, list)
            else x
        )
        listified = {k: listify(v) for k, v in wildcards.items()}
        missing_outputs = []
        complete_jobs = []
        for wc in dictproduct(listified):
            output, _ = self.rule.expand_output(wc)
            if self.checkpoints.created_output:
                missing_output = set(output) - set(self.checkpoints.created_output)
                if not missing_output:
                    complete_jobs.append(CheckpointJob(self.rule, output))
                else:
                    logger.debug(
                        f"Missing checkpoint output for {self.rule.name} "
                        f"(wildcards: {wc}): {','.join(missing_output)} of {','.join(output)}"
                    )
            missing_outputs.append(checkpoint_target(output[0]))
        if not missing_outputs:
            return complete_jobs
        raise IncompleteCheckpointException(self.rule, missing_outputs)


class CheckpointJob:
    __slots__ = ["rule", "output"]

    def __init__(self, rule: "Rule", output):
        self.output = output
        self.rule = rule
