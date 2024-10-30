from __future__ import annotations

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.rules import Rule


class OutputIndex:
    """Look up structure for rules, that can be queried by the output products which they create."""

    def __init__(self, rules: list[Rule]) -> None:
        self._entries = [
            (rule, str(o.constant_prefix()), str(o.constant_suffix()))
            for rule in rules
            for o in rule.products()
        ]

    def match(self, targetfile: str) -> set[Rule]:
        return {
            rule
            for rule, prefix, suffix in self._entries
            if targetfile.startswith(prefix) and targetfile.endswith(suffix)
        }

    def match_producer(self, targetfile: str) -> set[Rule]:
        return {rule for rule in self.match(targetfile) if rule.is_producer(targetfile)}
