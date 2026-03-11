__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2022, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
from abc import ABCMeta, abstractmethod
from pathlib import Path

from snakemake.caching.hash import ProvenanceHashMap
from snakemake.exceptions import (
    CacheMissException,
    MissingOutputFileCachePathException,
    WorkflowError,
)
from snakemake.io import apply_wildcards
from snakemake.jobs import Job
from snakemake.logging import logger

LOCATION_ENVVAR = "SNAKEMAKE_OUTPUT_CACHE"


class AbstractOutputFileCache:
    __metaclass__ = ABCMeta

    def __init__(self):
        try:
            self.cache_location = os.environ[LOCATION_ENVVAR]
        except KeyError:
            raise MissingOutputFileCachePathException()
        self.provenance_hash_map = ProvenanceHashMap()

    @abstractmethod
    async def store(self, job: Job, cache_mode):
        pass

    @abstractmethod
    async def fetch(self, job: Job, cache_mode):
        pass

    @abstractmethod
    async def exists(self, job: Job):
        pass

    def get_outputfiles(self, job: Job):
        if job.rule.output[0].is_multiext:
            prefix_len = len(
                apply_wildcards(job.rule.output[0].multiext_prefix, job.wildcards)
            )
            yield from ((f, f[prefix_len:]) for f in job.output)
        else:
            assert len(job.output) == 1, (
                "bug: multiple output files in cacheable job but multiext not used for declaring them"
            )
            # It is crucial to distinguish cacheable objects by the file extension.
            # Otherwise, for rules that generate different output based on the provided
            # extension a wrong cache entry can be returned.
            # Another nice side effect is that the cached files become more accessible
            # because their extension is presented in the cache dir.
            ext = Path(job.output[0]).suffix
            yield (job.output[0], ext)

    def raise_write_error(self, entry, exception=None):
        raise WorkflowError(
            f"Given output cache entry {entry} ($SNAKEMAKE_OUTPUT_CACHE={self.cache_location}) is not writeable.",
            *[exception],
        )

    def raise_read_error(self, entry, exception=None):
        raise WorkflowError(
            f"Given output cache entry {entry} ($SNAKEMAKE_OUTPUT_CACHE={self.cache_location}) is not readable.",
            *[exception],
        )

    def raise_cache_miss_exception(self, job):
        raise CacheMissException(f"Job {job} not yet cached.")
