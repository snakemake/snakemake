__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2019, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from snakemake.caching.hash import ProvenanceHashMap
from snakemake.caching import AbstractOutputFileCache
from snakemake.jobs import Job
from snakemake.io import get_flag_value, IOFile


class OutputFileCache(AbstractOutputFileCache):
    """
    A cache for output files that uses a provenance hash value that
    describes all steps, parameters, and software needed to generate
    each output file. This is the remote version.
    """

    def __init__(self, remote_provider):
        super().__init__()
        self.remote_provider = remote_provider

    def store(self, job: Job):
        outputfile = self.get_outputfile(job)
        entry = self._get_remote(job)

        # upload to remote
        entry.upload()

    def fetch(self, job: Job):
        self.check_job(job)
        entry = self._get_remote(job)
        outputfile = job.output[0]

        if not entry.exists():
            self.raise_cache_miss_exception(job)

        # download to outputfile
        entry.download()

    def exists(self, job: Job):
        self.check_job(job)
        outputfile = job.output[0]
        entry = self._get_remote(job)
        entry.exists()

    def _get_remote(self, job: Job):
        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        f = self.remote_provider.remote(
            "{}/{}".format(self.cache_location, provenance_hash)
        )
        remote = get_flag_value(f, "remote_object")
        remote._iofile = job.output[0]
        return remote
