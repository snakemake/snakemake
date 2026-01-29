__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2022, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import ABCMeta, abstractmethod
import os
from threading import Lock
from urllib.parse import quote

from snakemake.jobs import Job
from snakemake.io import apply_wildcards
from snakemake.exceptions import (
    MissingOutputFileCachePathException,
    WorkflowError,
    CacheMissException,
)
from snakemake.caching.hash import ProvenanceHashMap

LOCATION_ENVVAR = "SNAKEMAKE_OUTPUT_CACHE"


class AbstractOutputFileCache:
    __metaclass__ = ABCMeta

    def __init__(self):
        try:
            self.cache_location = os.environ[LOCATION_ENVVAR]
        except KeyError:
            raise MissingOutputFileCachePathException()
        self.provenance_hash_map = ProvenanceHashMap()

        # Locking mechanism to prevent concurrent execution of same cache key
        # Set of provenance hashes currently being computed
        self._active_hashes = set()
        # Single lock to protect the set
        self._active_hashes_lock = Lock()

    @abstractmethod
    async def store(self, job: Job, cache_mode):
        pass

    @abstractmethod
    async def fetch(self, job: Job, cache_mode):
        pass

    @abstractmethod
    async def exists(self, job: Job):
        pass

    def mark_if_schedulable(self, job: Job, cache_mode) -> bool:
        """Make a job as scheduled unless its provenance hash is already active

        Args:
            job (Job): Job to try scheduling
            cache_mode (str): Cache mode for this job

        Returns:
            bool: True if the job can be scheduled (hash was not marked yet).
                  False if the hash is already being computed by another job.
        """
        provenance_hash = self.provenance_hash_map.get_provenance_hash(job, cache_mode)
        with self._active_hashes_lock:
            if provenance_hash in self._active_hashes:
                # Hash already being computed by another job
                return False
            # Hash available, mark as being computed
            self._active_hashes.add(provenance_hash)
            return True

    def discard_mark(self, job: Job, cache_mode):
        """Release mark for a job.

        Called when job not selected by job_selector or after completion.

        Args:
            job (Job): Job whose mark should be released
            cache_mode (str): Cache mode for this job
        """
        provenance_hash = self.provenance_hash_map.get_provenance_hash(job, cache_mode)
        with self._active_hashes_lock:
            self._active_hashes.discard(provenance_hash)

    def get_outputfiles(self, job: Job):
        for _, fn in job.rule.output._allitems():
            yield apply_wildcards(fn, job.wildcards), quote(fn, safe="{}")

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

    def raise_cache_miss_exception(self, job):
        raise CacheMissException(f"Job {job} not yet cached.")
