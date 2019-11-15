__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2019, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from tempfile import TemporaryDirectory
from pathlib import Path
import os
import shutil
import stat

from snakemake.logging import logger
from snakemake.jobs import Job
from snakemake.exceptions import WorkflowError
from snakemake.caching.hash import ProvenanceHashMap
from snakemake.caching import LOCATION_ENVVAR, AbstractOutputFileCache


class OutputFileCache(AbstractOutputFileCache):
    """
    A cache for output files that uses a provenance hash value that
    describes all steps, parameters, and software needed to generate
    each output file.
    """

    def __init__(self):
        super().__init__()
        self.path = Path(self.cache_location)

    def check_writeable(self, entry):
        if not (os.access(self.path, os.W_OK) or os.access(self.path / entry, os.W_OK)):
            self.raise_write_error(entry)

    def check_readable(self, entry):
        if not os.access(self.path / entry, os.R_OK):
            self.raise_read_error(entry)

    def store(self, job: Job):
        """
        Store generated job output in the cache.
        """
        outputfile = self.get_outputfile(job, check_exists=True)

        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        self.check_writeable(provenance_hash)
        path = self.path / provenance_hash

        logger.info("Moving output file {} to cache.".format(outputfile))
        with TemporaryDirectory(dir=self.path) as tmpdirname:
            tmp = Path(tmpdirname) / provenance_hash
            # First move is performed into a tempdir (it might involve a copy if not on the same FS).
            # This is important, such that network filesystem latency
            # does not lead to concurrent writes to the same file.
            # We can use the plain copy method of shutil, because we do not care about the metadata.
            shutil.move(outputfile, tmp, copy_function=shutil.copy)
            # make readable/writeable for all
            os.chmod(
                tmp,
                stat.S_IRUSR
                | stat.S_IWUSR
                | stat.S_IRGRP
                | stat.S_IWGRP
                | stat.S_IROTH
                | stat.S_IWOTH,
            )

            # Move to the actual path (now we are on the same FS, hence move is atomic).
            # Here we use the default copy function, also copying metadata (which is important here).
            # It will always work, because we are guaranteed to be in the same FS.
            shutil.move(tmp, path)
        # now restore the outputfile via a symlink
        self.symlink(path, outputfile, utime=False)

    def fetch(self, job: Job):
        """
        Retrieve cached output file and copy to the place where the job expects it's output.
        """
        self.check_job(job)

        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        path = self.path / provenance_hash

        if not path.exists():
            self.raise_cache_miss_exception(job)

        self.check_readable(provenance_hash)

        outputfile = job.output[0]

        self.symlink(path, outputfile)

    def exists(self, job: Job):
        """
        Return True if job is already cached
        """
        self.check_job(job)

        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        path = self.path / provenance_hash

        if not path.exists():
            return False

        self.check_readable(provenance_hash)
        return True

    def symlink(self, path, outputfile, utime=True):
        if os.utime in os.supports_follow_symlinks or not utime:
            logger.info("Symlinking output file {} from cache.".format(outputfile))
            os.symlink(path, outputfile)
            if utime:
                os.utime(outputfile, follow_symlinks=False)
        else:
            logger.info(
                "Copying output file {} from cache (OS does not support updating the modification date of symlinks).".format(
                    outputfile
                )
            )
            shutil.copyfile(path, outputfile)
