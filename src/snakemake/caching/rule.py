# keep until 3.14
from __future__ import annotations

from dataclasses import dataclass
from enum import Flag, auto
from typing import TYPE_CHECKING

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger

if TYPE_CHECKING:
    from typing import Any, Optional, Self

    from snakemake.jobs import Job
    from snakemake.rules import Rule


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

        if len(rule.output) == 0:
            raise WorkflowError(
                "Rules without output files cannot be cached.", rule=rule
            )
        if len(rule.output) > 1:
            if not all(out.is_multiext for out in rule.output):
                raise WorkflowError(
                    "Rule is marked for between workflow caching but has multiple output files. "
                    "This is only allowed if multiext() is used to declare them (see docs on between "
                    "workflow caching).",
                    rule=rule,
                )
        if not rule.workflow.workflow_settings.cache:
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

    async def exists(self, job: Job):
        assert self.rule.workflow.output_file_cache is not None
        return await self.rule.workflow.output_file_cache.exists(job)

    async def store(self, job: Job):
        assert self.rule.workflow.output_file_cache is not None
        return await self.rule.workflow.output_file_cache.store(job)

    async def fetch(self, job: Job):
        assert self.rule.workflow.output_file_cache is not None
        return await self.rule.workflow.output_file_cache.fetch(job)
