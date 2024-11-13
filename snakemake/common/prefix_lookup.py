import bisect
from collections import defaultdict
from typing import Sequence, TypeVar, Generator

V = TypeVar("V")


class PrefixLookup:
    """A data structure for efficiently looking up entries by their prefix."""

    def __init__(self, entries: Sequence[tuple[str, V]]) -> None:
        grouped = defaultdict(list)
        for key, value in entries:
            grouped[key].append(value)
        self._entries = sorted(grouped.items(), key=lambda x: x[0])

    def match(self, key: str) -> set[V]:
        """Returns a set of all entries which are prefixes of the given key."""
        return set(self.match_iter(key))

    def match_iter(self, key: str) -> Generator[V, None, None]:
        """Yields all entries which are prefixes of the given key.

        E.g. if "abc" is the key then "abc", "ab", "a", and "" are considered valid prefixes.
        """
        stop_idx = bisect.bisect_right(self._entries, key, key=lambda x: x[0])
        last_len = float("inf")
        for index in range(stop_idx - 1, -1, -1):
            k, entries = self._entries[index]
            if key.startswith(k):
                for entry in entries:
                    yield entry
                last_len = len(k)
            elif len(k) > last_len:
                # If we find a non-matching prefix shorter than our last match, we can stop
                break
