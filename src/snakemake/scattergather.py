from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field
from math import floor
from typing import Any

from snakemake.exceptions import WorkflowError


@dataclass(frozen=True)
class ScatterProcess:
    items: int
    tolerance: float = 0.0
    impatient: bool = False

    def __post_init__(self):
        if not isinstance(self.items, int) or self.items <= 0:
            raise WorkflowError(
                "Scatter process item count has to be a positive integer."
            )
        if not isinstance(self.tolerance, (int, float)):
            raise WorkflowError("Scatter process tolerance has to be numeric.")
        if not 0.0 <= float(self.tolerance) <= 1.0:
            raise WorkflowError(
                "Scatter process tolerance has to be between 0.0 and 1.0."
            )
        object.__setattr__(self, "tolerance", float(self.tolerance))
        if not isinstance(self.impatient, bool):
            raise WorkflowError("Scatter process impatient flag has to be boolean.")

    @property
    def max_failures(self) -> int:
        # Tolerance is expressed as a fraction of the total number of scatter items.
        return floor(self.items * self.tolerance)

    @property
    def required_successes(self) -> int:
        # Impatient mode can release gather once this many branches have succeeded.
        return self.items - self.max_failures

    def with_items(self, items: int) -> ScatterProcess:
        return ScatterProcess(
            items=items, tolerance=self.tolerance, impatient=self.impatient
        )

    @classmethod
    def parse(cls, value: Any) -> ScatterProcess:
        # Keep the legacy integer form while accepting richer process definitions.
        if isinstance(value, cls):
            return value
        if isinstance(value, int):
            return cls(items=value)
        if isinstance(value, Mapping):
            supported_keys = {"items", "tolerance", "impatient"}
            unknown = set(value) - supported_keys
            if unknown:
                raise WorkflowError(
                    "Unsupported scatter process options: "
                    + ", ".join(sorted(map(str, unknown)))
                )
            if "items" not in value:
                raise WorkflowError(
                    "Scatter process mappings require an 'items' entry."
                )
            if "impatient" in value and "tolerance" not in value:
                raise WorkflowError(
                    "Scatter process 'impatient' mode requires a 'tolerance' bound."
                )
            return cls(
                items=value["items"],
                tolerance=value.get("tolerance", 0.0),
                impatient=value.get("impatient", False),
            )
        raise WorkflowError(
            "Scatter process definitions have to be integers, mappings, or "
            "scattergather_process(...) objects."
        )


def scattergather_process(
    items: int, tolerance: float = 0.0, impatient: bool = False
) -> ScatterProcess:
    return ScatterProcess(items=items, tolerance=tolerance, impatient=impatient)


@dataclass(frozen=True)
class ScatterItem:
    process: str
    item: int
    total: int
    role: str


@dataclass
class ScatterProcessState:
    definition: ScatterProcess
    gather_jobs: set[Any] = field(default_factory=set)
    # Upstream jobs that can decide the fate of a given scatter item.
    producers: dict[int, set[Any]] = field(default_factory=lambda: defaultdict(set))
    successes: set[int] = field(default_factory=set)
    failures: set[int] = field(default_factory=set)
    ignored: set[int] = field(default_factory=set)
    # In impatient mode, freeze the accepted branches once enough succeeded.
    accepted: frozenset[int] | None = None

    @property
    def all_items(self) -> set[int]:
        return set(range(1, self.definition.items + 1))

    @property
    def active_items(self) -> set[int]:
        """
        The set of scatter items that are still active
        (i.e. those which could still proceed to the gather stage).
        """
        if self.accepted is not None:
            # In 'impatient' mode, only keep the explicitly accepted branches.
            return set(self.accepted)
        return self.all_items - self.failures


@dataclass(frozen=True)
class ScatterFailureResult:
    tolerated_processes: tuple[str, ...] = ()
    fatal_processes: tuple[str, ...] = ()
    ignored_processes: tuple[str, ...] = ()

    @property
    def handled(self) -> bool:
        return bool(
            self.tolerated_processes or self.fatal_processes or self.ignored_processes
        )

    @property
    def is_fatal(self) -> bool:
        return bool(self.fatal_processes)
