__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2019, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import ABCMeta, abstractmethod
import os

from snakemake.jobs import Job
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

    def get_outputfile(self, job: Job):
        output = list(
            job.expanded_output
        )  # TODO remove one dynamic is removed from codebase
        assert len(output) == 1, "Bug: Only single output files are supported."
        outputfile = output[0]
        assert os.path.exists(
            outputfile
        ), "Bug: Output file does not exist although it should be cached."
        return outputfile

    def check_job(self, job: Job):
        assert len(job.output) == 1, "Bug: Only single output files are supported."

    def raise_write_error(self, entry):
        raise WorkflowError(
            "Given output cache entry {} ($SNAKEMAKE_OUTPUT_CACHE={}) is not writeable.".format(
                entry, self.cache_location
            )
        )

    def raise_read_error(self, entry):
        raise WorkflowError(
            "Given output cache entry {} ($SNAKEMAKE_OUTPUT_CACHE={}) is not readable.".format(
                entry, self.cache_location
            )
        )

    def raise_cache_miss_exception(self, job):
        raise CacheMissException("Job {} not yet cached.".format(job))
