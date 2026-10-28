# keep until 3.14
from __future__ import annotations

import functools
from dataclasses import dataclass
from enum import Flag, auto
from typing import TYPE_CHECKING

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger

if TYPE_CHECKING:
    from typing import Any, Optional, Self

    from snakemake.jobs import Job
    from snakemake.rules import Rule


def _requires_cache(method):
    @functools.wraps(method)
    async def wrapper(self, *args, **kwargs):
        if self.rule.workflow.output_file_cache is None:
            raise RuntimeError(
                f"Cache backend unavailable for rule '{self.rule}' during"
                f" '{method.__name__}' operation."
            )
        return await method(self, *args, **kwargs)

    return wrapper


class CacheFlag(Flag):
    output = auto()
    omit_software = auto()
    omit_storage_content = auto()


@dataclass
class RuleCache:
    rule: Rule
    flag: CacheFlag

    @classmethod
    def from_rule(cls, rule: Rule, cache: Any) -> Optional[Self]:
        if not cache:
            flag = CacheFlag(0)
        elif cache is True or cache == "all":
            flag = CacheFlag.output
        else:
            flag = CacheFlag(0)
            for x in [cache] if isinstance(cache, str) else cache:
                if x == "omit-software":
                    flag |= CacheFlag.output | CacheFlag.omit_software
                elif x == "hash-omit-software":
                    flag |= CacheFlag.omit_software
                elif x == "omit-storage-content":
                    flag |= CacheFlag.output | CacheFlag.omit_storage_content
                elif x == "hash-omit-storage-content":
                    flag |= CacheFlag.omit_storage_content
                else:
                    raise WorkflowError(
                        "Invalid value for cache directive. Use 'all', "
                        "'omit-software', 'hash-omit-software', 'omit-storage-content', "
                        "'hash-omit-storage-content'.",
                        rule=rule,
                    )

        return cls(rule, flag)

    @property
    def output(self):
        return CacheFlag.output in self.flag

    @property
    def omit_software(self):
        return CacheFlag.omit_software in self.flag

    @property
    def omit_storage_content(self):
        return CacheFlag.omit_storage_content in self.flag

    def check(self):
        if not self.output:
            return

        rule = self.rule

        if not rule.output:
            raise WorkflowError(
                "Rules without output files cannot be cached.", rule=rule
            )
        if len(rule.output) > 1:
            if not all(out.is_multiext for out in rule.output) and not all(
                name is not None for name, _ in rule.output._allitems()
            ):
                raise WorkflowError(
                    "Rule is marked for between workflow caching but has multiple output files. "
                    "This is only allowed if multiext() is used to declare them, or if all "
                    "output files are named (see docs on between workflow caching).",
                    rule=rule,
                )
        if rule.workflow.workflow_settings.cache is None:
            logger.warning(
                f"Workflow defines that rule {rule.name} is eligible for caching between workflows "
                "(use the --cache argument to enable this)."
            )

        if rule.benchmark:
            raise WorkflowError(
                "Rules with a benchmark directive may not be marked as eligible "
                "for between-workflow caching at the same time. The reason is that "
                "when the result is taken from cache, there is no way to fill the benchmark file with "
                "any reasonable values. Either remove the benchmark directive or disable "
                "between-workflow caching for this rule.",
                rule=rule,
            )

    @_requires_cache
    async def exists(self, job: Job):
        return await self.rule.workflow.output_file_cache.exists(job)

    @_requires_cache
    async def store(self, job: Job):
        return await self.rule.workflow.output_file_cache.store(job)

    @_requires_cache
    async def fetch(self, job: Job):
        return await self.rule.workflow.output_file_cache.fetch(job)
