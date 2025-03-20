from enum import Enum
from typing import Callable, Iterable, Union
from snakemake.io import flag

VALID_PATTERNS = frozenset(("random", "sequential", "sequential-multi"))


class AccessPattern(Enum):
    RANDOM = "random"
    SEQUENTIAL = "sequential"
    MULTI = "multi"

    def __str__(self):
        return self.value


item_type = Union[str, Iterable[str], Callable, Iterable[Callable]]


class AccessPatternFactory:
    @classmethod
    def random(cls, item: item_type):
        return flag(item, "access_pattern", AccessPattern.RANDOM)

    @classmethod
    def sequential(cls, item: item_type):
        return flag(item, "access_pattern", AccessPattern.SEQUENTIAL)

    @classmethod
    def multi(cls, item: item_type):
        return flag(item, "access_pattern", AccessPattern.MULTI)
