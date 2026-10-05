from typing import TYPE_CHECKING
from snakemake.exceptions import IncompleteCheckpointException, WorkflowError
from snakemake.io import _IOFile, checkpoint_target, get_flag_value
from snakemake.logging import logger

if TYPE_CHECKING:
    from snakemake.rules import Rule
    from snakemake.iocontainers import OutputFiles
    from snakemake.io.typed import _TypedFile


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


class CheckpointTypedFile(_IOFile):
    """A resolved checkpoint output with synchronous typed content access."""

    def load(self):
        """Synchronously read this checkpoint output using its typed loader."""
        x = self.plainstr
        if TYPE_CHECKING:
            assert isinstance(x, _TypedFile)
        return x.load()


def typed_guard(iofile: _IOFile):
    typed_factory = get_flag_value(iofile, "typed")
    if typed_factory:
        return CheckpointTypedFile(iofile, iofile.rule)
    return iofile


class CheckpointJob:
    __slots__ = ["rule", "output"]

    def __init__(self, rule: "Rule", output: "OutputFiles"):
        # Keep file flags and names when outputs are reused as downstream inputs.
        self.output = output.__class__(toclone=output, custom_map=typed_guard)
        self.rule = rule
