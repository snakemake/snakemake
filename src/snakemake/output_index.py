from __future__ import annotations

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

from typing import TYPE_CHECKING

from snakemake.common.prefix_lookup import PrefixLookup

if TYPE_CHECKING:
    from snakemake.rules import Rule


class OutputIndex:
    """Look up structure for rules, that can be queried by the output products which they create."""

    def __init__(self, rules: list[Rule]) -> None:
        entries = []
        for rule in rules:
            for product in rule.products():
                prefix = str(product.constant_prefix())
                suffix = str(product.constant_suffix())
                entries.append((prefix, (rule, suffix)))
        self._lookup = PrefixLookup(entries=entries)

    def match(self, targetfile: str) -> set[Rule]:
        """Returns all rules that match the given target file, considering only the prefix and suffix up to the
        first wildcard.

        To further verify the match, the returned rules should be checked with ``Rule.is_producer(targetfile)``.
        """
        return {
            rule
            for prefix, (rule, suffix) in self._lookup.match_iter(targetfile)
            if targetfile.endswith(suffix)
        }

    def match_producers(self, targetfile: str) -> set[Rule]:
        """Returns all rules that match and produce the given target file."""
        return {rule for rule in self.match(targetfile) if rule.is_producer(targetfile)}
