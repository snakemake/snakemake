__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2021, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import ABCMeta, abstractmethod
import os

from snakemake.jobs import Job
from snakemake.io import is_flagged, get_flag_value, apply_wildcards
from snakemake.exceptions import WorkflowError, CacheMissException
from snakemake.caching.hash import ProvenanceHashMap

LOCATION_ENVVAR = "SNAKEMAKE_OUTPUT_CACHE"


class AbstractOutputFileCache:
    __metaclass__ = ABCMeta

    def __init__(self):
        try:
            self.cache_location = os.environ[LOCATION_ENVVAR]
        except KeyError:
            raise WorkflowError(
                "Output file cache activated (--cache), but no cache "
                "location specified. Please set the environment variable "
                "${}.".format(LOCATION_ENVVAR)
            )
        self.provenance_hash_map = ProvenanceHashMap()

    @abstractmethod
    def store(self, job: Job):
        pass

    @abstractmethod
    def fetch(self, job: Job):
        pass

    @abstractmethod
    def exists(self, job: Job):
        pass

    def get_outputfiles(self, job: Job):
        if job.rule.output[0].is_multiext:
            prefix_len = len(
                apply_wildcards(job.rule.output[0].multiext_prefix, job.wildcards)
            )
            yield from ((f, f[prefix_len:]) for f in job.output)
        else:
            yield (job.output[0], "")

    def raise_write_error(self, entry, exception=None):
        raise WorkflowError(
            "Given output cache entry {} ($SNAKEMAKE_OUTPUT_CACHE={}) is not writeable.".format(
                entry, self.cache_location
            ),
            *[exception],
        )

    def raise_read_error(self, entry, exception=None):
        raise WorkflowError(
            "Given output cache entry {} ($SNAKEMAKE_OUTPUT_CACHE={}) is not readable.".format(
                entry, self.cache_location
            ),
            *[exception],
        )

    def raise_cache_miss_exception(self, job):
        raise CacheMissException("Job {} not yet cached.".format(job))
