from contextvars import ContextVar, Token
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

    def __getattr__(self, name):
        if name in self._rules:
            return self._rules[name]
        if name == "waitfor":
            # Alias for `with checkpoints.waitfor as waitfor:`
            # to make this semantically easier to read
            return self
        raise WorkflowError(
            f"Checkpoint {name} is not defined in this workflow. "
            f"Available checkpoints: {', '.join(self._rules)}"
        )

    def register(self, rule: "Rule", fallback_name=None):
        self._rules[fallback_name or rule.name] = Checkpoint(rule, self._parent)

    def __enter__(self):
        """Return a :class:`CheckpointsWaitforProxy` `with`in which multiple
          checkpoint dependencies can be declared at once
          instead of raising :class:`IncompleteCheckpointException` immediately.

        The returned proxy is registered as the active collector via a module-level
        :class:`~contextvars.ContextVar`, which isolates state per-coroutine and
        per-thread without modifying any shared state on :class:`CheckpointsProxy` itself.
        The token from :meth:`~contextvars.ContextVar.set` is stored on the proxy so
        :meth:`__exit__` can restore the previous collector, correctly supporting
        nested `with checkpoints as waitfor:` blocks.
        """
        waitfor = CheckpointsWaitforProxy(self)
        waitfor._token = _current_collector.set(waitfor)
        return waitfor

    def __exit__(self, exc_type, exc, tb):
        """
        Restore the previous collector and raise collected checkpoint exceptions.

        Resets the :class:`~contextvars.ContextVar` to its value before the
        corresponding :meth:`__enter__`, ensuring correct behavior for nested
        `with checkpoints as waitfor:` blocks.

        Three cases on exit:
        - :class:`IncompleteCheckpointException` escaped
          (e.g. user called `checkpoints.foo` instead of `waitfor.foo` inside the block):
            folded back into the cache so all pending dependencies are reported together.
        - No exception: if any :class:`IncompleteCheckpointException` in cache,
          the first is raised with the rest attached as `.extra`.
        - Unrelated exception: left to propagate normally.
        """
        waitfor: CheckpointsWaitforProxy = _current_collector.get()
        _current_collector.reset(waitfor._token)
        cache = waitfor.cache
        if exc_type is IncompleteCheckpointException:
            # User accidentally used `checkpoints.foo` instead of `waitfor.foo`.
            # Fold the escaped exception into the cache so all dependencies are reported together.
            cache.append(exc)
            exc_type = None
        if exc_type is None and cache:
            e, *_e = cache
            e.extra = _e
            raise e
        return False


# It is safe to use a single module-level ContextVar for all CheckpointsProxy instances.
# 1. See: https://docs.python.org/3/library/contextvars.html:
#    > ContextVar should be created at the top module level, never in closures.
# 2. Although multiple `CheckpointsProxy` instances (one per *workflow module*) share this variable,
#    their `with` blocks are never interleaved within the same coroutine:
#    each block completes before the next checkpoint function is evaluated.
_current_collector: ContextVar = ContextVar("checkpoint_proxy_waitfor", default=None)


class CheckpointsWaitforProxy:
    """Collect :class:`IncompleteCheckpointException` into :attr:`cache`
    instead of raising them immediately,
    so multiple checkpoint dependencies can be declared at once during DAG construction.
    """

    def __init__(self, parent: CheckpointsProxy):
        self.parent = parent
        self.cache: list[IncompleteCheckpointException] = []
        self._token: Token

    def __getattr__(self, name):
        checkpoint = getattr(self.parent, name)
        return CheckpointLater(checkpoint, self.cache)


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
