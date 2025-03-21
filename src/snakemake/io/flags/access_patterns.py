from enum import Enum
from snakemake.io import flag
from snakemake.io.flags import FlaggableItemOrIterable


class AccessPattern(Enum):
    RANDOM = "random"
    SEQUENTIAL = "sequential"
    MULTI = "multi"

    def __str__(self):
        return self.value


class AccessPatternFactory:
    @classmethod
    def random(cls, item: FlaggableItemOrIterable):
        return flag(item, "access_pattern", AccessPattern.RANDOM)

    @classmethod
    def sequential(cls, item: FlaggableItemOrIterable):
        return flag(item, "access_pattern", AccessPattern.SEQUENTIAL)

    @classmethod
    def multi(cls, item: FlaggableItemOrIterable):
        return flag(item, "access_pattern", AccessPattern.MULTI)
