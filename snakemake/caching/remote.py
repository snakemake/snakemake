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
        entry = self._get_remote(job, check_output_exists=True)

        # upload to remote
        try:
            entry.upload()
        except Exception as e:
            self.raise_write_error(entry, exception=e)

    def fetch(self, job: Job):
        self.check_job(job)
        entry = self._get_remote(job)

        if not entry.exists():
            self.raise_cache_miss_exception(job)

        # download to outputfile
        try:
            entry.download()
        except Exception as e:
            self.raise_read_error(entry, exception=e)

    def exists(self, job: Job):
        self.check_job(job)
        entry = self._get_remote(job)

        try:
            return entry.exists()
        except Exception as e:
            self.raise_read_error(entry, exception=e)

    def _get_remote(self, job: Job, check_output_exists=False):
        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        f = self.remote_provider.remote(
            "{}/{}".format(self.cache_location, provenance_hash)
        )
        remote = get_flag_value(f, "remote_object")
        remote._iofile = self.get_outputfile(job, check_exists=check_output_exists)
        return remote
