__authors__ = "Johannes Köster, Sven Nahnsen"
__copyright__ = "Copyright 2021, Johannes Köster, Sven Nahnsen"
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
        # make readable/writeable for all
        self.file_permissions = (
            stat.S_IRUSR
            | stat.S_IWUSR
            | stat.S_IRGRP
            | stat.S_IWGRP
            | stat.S_IROTH
            | stat.S_IWOTH
        )
        # directories need to have exec permission as well (for opening)
        self.dir_permissions = (
            self.file_permissions | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        )

    def check_writeable(self, cachefile):
        if not (os.access(cachefile.parent, os.W_OK) or os.access(cachefile, os.W_OK)):
            self.raise_write_error(cachefile)

    def check_readable(self, cachefile):
        if not os.access(cachefile, os.R_OK):
            self.raise_read_error(cachefile)

    def store(self, job: Job):
        """
        Store generated job output in the cache.
        """

        if not os.access(self.path, os.W_OK):
            raise WorkflowError(
                "Cannot access cache location {}. Please ensure that "
                "it is present and writeable.".format(self.path)
            )
        with TemporaryDirectory(dir=self.path) as tmpdirname:
            tmpdir = Path(tmpdirname)

            for outputfile, cachefile in self.get_outputfiles_and_cachefiles(job):
                if not os.path.exists(outputfile):
                    raise WorkflowError(
                        "Cannot move output file {} to cache. It does not exist "
                        "(maybe it was not created by the job?)."
                    )
                self.check_writeable(cachefile)
                logger.info("Moving output file {} to cache.".format(outputfile))

                tmp = tmpdir / cachefile.name
                # First move is performed into a tempdir (it might involve a copy if not on the same FS).
                # This is important, such that network filesystem latency
                # does not lead to concurrent writes to the same file.
                # We can use the plain copy method of shutil, because we do not care about the metadata.
                shutil.move(outputfile, tmp, copy_function=shutil.copy)

                self.set_permissions(tmp)

                # Move to the actual path (now we are on the same FS, hence move is atomic).
                # Here we use the default copy function, also copying metadata (which is important here).
                # It will always work, because we are guaranteed to be in the same FS.
                shutil.move(tmp, cachefile)
                # now restore the outputfile via a symlink
                self.symlink(cachefile, outputfile, utime=False)

    def fetch(self, job: Job):
        """
        Retrieve cached output file and symlink to the place where the job expects it's output.
        """
        for outputfile, cachefile in self.get_outputfiles_and_cachefiles(job):

            if not cachefile.exists():
                self.raise_cache_miss_exception(job)

            self.check_readable(cachefile)
            if cachefile.is_dir():
                # For directories, create a new one and symlink each entry.
                # Then, the .snakemake_timestamp of the new dir is touched
                # by the executor.
                outputfile.mkdir(parents=True, exist_ok=True)
                for f in cachefile.iterdir():
                    self.symlink(f, outputfile / f.name)
            else:
                self.symlink(cachefile, outputfile)

    def exists(self, job: Job):
        """
        Return True if job is already cached
        """
        for outputfile, cachefile in self.get_outputfiles_and_cachefiles(job):

            if not cachefile.exists():
                return False

            self.check_readable(cachefile)
        return True

    def get_outputfiles_and_cachefiles(self, job: Job):
        provenance_hash = self.provenance_hash_map.get_provenance_hash(job)
        base_path = self.path / provenance_hash

        return (
            (Path(outputfile), base_path.with_suffix(ext))
            for outputfile, ext in self.get_outputfiles(job)
        )

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

    def set_permissions(self, entry):
        # make readable/writeable for all
        if entry.is_dir():
            # recursively apply permissions for all contained files
            for root, dirs, files in os.walk(entry):
                root = Path(root)
                for d in dirs:
                    os.chmod(root / d, self.dir_permissions)
                for f in files:
                    os.chmod(root / f, self.file_permissions)
            os.chmod(entry, self.dir_permissions)
        else:
            os.chmod(entry, self.file_permissions)
