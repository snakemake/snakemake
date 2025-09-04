from dataclasses import dataclass, field
from typing import Callable, Iterable, Mapping, Set, Union

from snakemake.io import get_flag_store_keys, is_flagged


# TODO move all flag definitions to this module

FlaggableItem = Union[str, Callable]

FlaggableItemOrIterable = Union[FlaggableItem, Iterable[FlaggableItem]]


CONTRADICTING_FLAGS = {
    "temp": {"protected"},
    "protected": {"temp"},
}


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

    def apply(self, item: FlaggableItem):
        for flag_ in self:
            store_keys = self.get_store_keys(flag_)
            if any(
                is_flagged(item, store_key)
                or any(
                    is_flagged(item, contradicting)
                    for contradicting in self.get_contradicting_flags(store_key)
                )
                for store_key in store_keys
            ):
                continue
            item = flag_(item)
        return item

    @staticmethod
    def get_contradicting_flags(store_key: str):
        return CONTRADICTING_FLAGS.get(store_key, set())
