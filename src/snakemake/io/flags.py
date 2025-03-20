from dataclasses import dataclass, field
from typing import Callable, Mapping, Set

from snakemake.io import get_flag_store_keys


# TODO move all flag definitions to this module


@dataclass
class DefaultFlags:
    flags_to_store_keys: Mapping[Callable, Set[str]] = field(default_factory=dict)

    def register_flags(self, *flags: Callable):
        self.flags_to_store_keys.clear()
        for flag in flags:
            keys = get_flag_store_keys(flag)
            self.flags_to_store_keys[flag] = keys

    def __iter__(self):
        return iter(self.flags_to_store_keys)

    def clear(self):
        self.flags_to_store_keys.clear()

    def get_store_keys(self, flag: Callable):
        return self.flags_to_store_keys[flag]
