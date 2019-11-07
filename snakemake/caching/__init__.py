__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2019, Johannes Köster, Sven Nahnsen"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from tempfile import NamedTemporaryFile
from pathlib import Path
import os
import shutil

from snakemake.logging import logger
from snakemake.jobs import Job
from snakemake.exceptions import WorkflowError, CacheMissException
from snakemake.caching.hash import ProvenanceHashMap

LOCATION_ENVVAR = "SNAKEMAKE_OUTPUT_CACHE"


class OutputFileCache:
    """
    A cache for output files that uses a provenance hash value that
    describes all steps, parameters, and software needed to generate
    each output file.
    """

    def __init__(self):
        try:
            self.path = Path(os.environ[LOCATION_ENVVAR])
        except KeyError:
            raise WorkflowError(
                "Output file cache activated (--cache), but no cache "
                "location specified. Please set the environment variable "
                "${}.".format(LOCATION_ENVVAR)
            )
        self.provenance_hash_map = ProvenanceHashMap()

    def check_writeable(self, entry):
        if not (os.access(self.path, os.W_OK) or os.access(self.path / entry, os.W_OK)):

            raise WorkflowError(
                "Given output cache entry {} ($SNAKEMAKE_OUTPUT_CACHE={}) is not writeable.".format(
                    entry, self.path
                )
            )

    def check_readable(self, entry):
        if not os.access(self.path / entry, os.R_OK):
            raise WorkflowError(
                "Given output cache entry {} ($SNAKEMAKE_OUTPUT_CACHE={}) is not readable.".format(
                    entry, self.path
                )
            )

    def store(self, job: Job):
        """
        Store generated job output in the cache.
        """
        output = list(job.expanded_output)
        assert (
            len(output) == 1
        ), "Bug: Only single output files are supported"
        outputfile = output[0]

        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        self.check_writeable(provenance_hash)
        path = self.path / provenance_hash
        # copy output file
        assert os.path.exists(outputfile)

        logger.info("Copying output file {} to cache.".format(outputfile))
        with NamedTemporaryFile(dir=self.path, delete=False) as tmp, open(outputfile, "rb") as out:
            # Copy is performed into a tempfile.
            # This is important, such that network filesystem latency 
            # does not lead to concurrent writes to the same file.
            shutil.copyfileobj(out, tmp)
        
        # rename (as atomic as possible) to the actual path
        os.rename(tmp.name, path)

    def fetch(self, job: Job):
        """
        Retrieve cached output file and copy to the place where the job expects it's output.
        """
        assert len(job.output) == 1, "Bug: Only single output files are supported"

        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        path = self.path / provenance_hash

        if not path.exists():
            raise CacheMissException("Job {} not yet cached.".format(job))

        self.check_readable(provenance_hash)

        outputfile = job.output[0]
        
        if os.utime in os.supports_follow_symlinks:
            logger.info("Symlinking output file {} from cache.".format(outputfile))
            os.symlink(path, outputfile)
            os.utime(outputfile, follow_symlinks=False)
        else:
            logger.info("Copying output file {} from cache (OS does not support updating the modification date of symlinks).".format(outputfile))
            shutil.copyfile(path, outputfile)

    
    def exists(self, job: Job):
        """
        Return True if job is already cached
        """
        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        path = self.path / provenance_hash

        if not path.exists():
            return False

        self.check_readable(provenance_hash)
        return True