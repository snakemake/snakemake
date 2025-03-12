__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2022, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
from pathlib import Path

from snakemake.caching import AbstractOutputFileCache
from snakemake.exceptions import WorkflowError
from snakemake.jobs import Job
from snakemake_interface_storage_plugins.storage_provider import StorageProviderBase


class OutputFileCache(AbstractOutputFileCache):
    """
    A cache for output files that uses a provenance hash value that
    describes all steps, parameters, and software needed to generate
    each output file. This is the remote version.
    """

    def __init__(self, storage_provider: StorageProviderBase):
        super().__init__()
        self.storage_provider = storage_provider

    async def store(self, job: Job, cache_mode: str):
        for entry in self._get_storage_objects(
            job, cache_mode, check_output_exists=True
        ):
            # upload to remote
            try:
                await entry.managed_store()
            except Exception as e:
                self.raise_write_error(entry, exception=e)

    async def fetch(self, job: Job, cache_mode: str):
        for entry in self._get_storage_objects(job, cache_mode):
            if not await entry.managed_exists():
                self.raise_cache_miss_exception(job)

            # download to outputfile
            try:
                await entry.managed_retrieve()
            except Exception as e:
                self.raise_read_error(entry, exception=e)

    async def exists(self, job: Job, cache_mode: str):
        for entry in self._get_storage_objects(job, cache_mode):
            try:
                return await entry.managed_exists()
            except Exception as e:
                self.raise_read_error(entry, exception=e)

    def _get_storage_objects(
        self, job: Job, cache_mode: str, check_output_exists=False
    ):
        provenance_hash = self.provenance_hash_map.get_provenance_hash(job, cache_mode)

        for outputfile, ext in self.get_outputfiles(job):
            if check_output_exists and not os.path.exists(outputfile):
                raise WorkflowError(
                    "Cannot move output file {} to cache. It does not exist "
                    "(maybe it was not created by the job?)."
                )

            storage_object = self.storage_provider.object(
                f"{self.cache_location}/{provenance_hash}{ext}"
            )
            storage_object.set_local_path(Path(outputfile))
            yield storage_object
