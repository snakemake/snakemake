from enum import Enum
from snakemake.io import flag
from snakemake.io.flags import FlaggableItemOrIterable

from snakemake_interface_common.io import AccessPattern

STORE_KEY = "access_pattern"


class AccessPatternFactory:
    @classmethod
    def random(cls, item: FlaggableItemOrIterable):
        return flag(item, STORE_KEY, AccessPattern.RANDOM)

    @classmethod
    def sequential(cls, item: FlaggableItemOrIterable):
        return flag(item, STORE_KEY, AccessPattern.SEQUENTIAL)

    @classmethod
    def multi(cls, item: FlaggableItemOrIterable):
        return flag(item, STORE_KEY, AccessPattern.MULTI)
