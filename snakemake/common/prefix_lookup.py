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

    def match(self, query: str) -> set[V]:
        """Returns a set of all entry values which are attached to prefixes of the given key."""
        return set(m[1] for m in self.match_iter(query))

    def match_iter(self, query: str) -> Generator[V, None, None]:
        """Yields all entries which are prefixes of the given key.

        E.g. if "abc" is the key then "abc", "ab", "a", and "" are considered
        valid prefixes.

        Yields entries as (key, value) tuples as supplied to __init__()
        """
        stop_idx = len(self._entries)
        while stop_idx:
            stop_idx = bisect.bisect_right(
                self._entries, query, hi=stop_idx, key=lambda x: x[0]
            )

            k, entries = self._entries[stop_idx - 1]
            if query.startswith(k):
                for entry in entries:
                    yield (k, entry)

            if not query or not k:
                # Exit loop, if iteration has reached "" query or key
                break

            query = self._common_prefix(k, query)

    def _common_prefix(self, s1: str, s2: str) -> str:
        """Utility function that gives the common prefix of two strings.

        Except, if the strings are identical, returns the string minus one letter.
        Behaviour is undefined if either string is empty.
        """
        for i, (l1, l2) in enumerate(zip(s1, s2)):
            if l1 != l2:
                break

        return s1[:i]
