# keep until 3.14 to avoid circular imports
from __future__ import annotations

__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2022, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import ABCMeta, abstractmethod
import os
from pathlib import Path
from typing import TYPE_CHECKING

from snakemake.io import apply_wildcards
from snakemake.exceptions import (
    MissingOutputFileCachePathException,
    WorkflowError,
    CacheMissException,
)
from snakemake.caching.hash import ProvenanceHashMap

if TYPE_CHECKING:
    from snakemake.jobs import Job

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
    async def store(self, job: Job):
        pass

    @abstractmethod
    async def fetch(self, job: Job):
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
        elif len(job.output) == 1:
            # It is crucial to distinguish cacheable objects by the file extension.
            # Otherwise, for rules that generate different output based on the provided
            # extension a wrong cache entry can be returned.
            # Another nice side effect is that the cached files become more accessible
            # because their extension is presented in the cache dir.
            yield (job.output[0], Path(job.output[0]).suffix)
        else:
            for name, outputfile in job.output._allitems():
                assert (
                    name is not None
                ), "bug: multiple unnamed output files in cacheable job"
                yield (outputfile, f"_{name}{Path(outputfile).suffix}")

    def raise_write_error(self, entry, exception=None):
        raise WorkflowError(
            "Given output cache entry {} ($SNAKEMAKE_OUTPUT_CACHE={}) is not writeable.".format(
                entry, self.cache_location
            ),
            exception,
        )

    def raise_read_error(self, entry, exception=None):
        raise WorkflowError(
            "Given output cache entry {} ($SNAKEMAKE_OUTPUT_CACHE={}) is not readable.".format(
                entry, self.cache_location
            ),
            exception,
        )

    def raise_cache_miss_exception(self, job: Job):
        raise CacheMissException(f"Job {job} not yet cached.")
