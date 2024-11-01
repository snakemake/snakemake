from __future__ import annotations

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

import bisect
from collections import defaultdict

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.rules import Rule


class OutputIndex:
    """Look up structure for rules, that can be queried by the output products which they create."""

    def __init__(self, rules: list[Rule]) -> None:
        entries = defaultdict(list)
        for rule in rules:
            for product in rule.products():
                prefix = str(product.constant_prefix())
                suffix = str(product.constant_suffix())
                entries[prefix].append((rule, suffix))
        self._entries = sorted(entries.items())

    def match(self, targetfile: str) -> set[Rule]:
        """Returns all rules that match the given target file, considering only the prefix and suffix up to the
        first wildcard.

        To further verify the match, the returned rules should be checked with ``Rule.is_producer(targetfile)``.
        """
        stop_idx = bisect.bisect_right(self._entries, targetfile, key=lambda x: x[0])
        hits = set()

        first = True
        for key, entries in self._entries[:stop_idx][::-1]:
            if targetfile.startswith(key):
                hits.update(rule for rule, suffix in entries if targetfile.endswith(suffix))
                first = False
            elif not first:
                break

        return hits

    def match_producer(self, targetfile: str) -> set[Rule]:
        return {rule for rule in self.match(targetfile) if rule.is_producer(targetfile)}
