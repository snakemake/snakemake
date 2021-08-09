__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import hashlib
from pathlib import Path
import re
import os
import tempfile
import io

from snakemake.common import is_local_file, get_appdirs
from snakemake.exceptions import WorkflowError
from snakemake.io import git_content


# TODO also use sourcecache for script and wrapper code!


class SourceCache:
    cache_whitelist = [
        "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/\d+\.\d+.\d+"
    ]  # TODO add more prefixes for uris that are save to be cached

    def __init__(self):
        self.cache = Path(
            os.path.join(get_appdirs().user_cache_dir, "snakemake/source-cache")
        )
        self.runtime_cache = tempfile.TemporaryDirectory(
            suffix="snakemake-runtime-source-cache"
        )
        self.cacheable_prefixes = re.compile("|".join(self.cache_whitelist))

    def lock_cache(self, entry):
        from filelock import FileLock

        return FileLock(entry.with_suffix(".lock"))

    def is_persistently_cacheable(self, path_or_uri):
        # TODO remove special git url handling once included in smart_open
        if path_or_uri.startswith("git+file:"):
            return False
        return is_local_file(path_or_uri) and self.cacheable_prefixes.match(path_or_uri)

    def open(self, path_or_uri, mode="r"):
        cache_entry = self._cache(path_or_uri)
        return self._open(cache_entry, mode)

    def get_path(self, path_or_uri, mode="r"):
        cache_entry = self._cache(path_or_uri)
        return cache_entry

    def _cache_entry(self, path_or_uri):
        urihash = hashlib.sha256()
        urihash.update(path_or_uri.encode())
        urihash = urihash.hexdigest()

        # TODO add git support to smart_open!
        if self.is_persistently_cacheable(path_or_uri):
            # check cache
            return self.cache / urihash
        else:
            # check runtime cache
            return Path(self.runtime_cache.name) / urihash

    def _cache(self, path_or_uri):
        cache_entry = self._cache_entry(path_or_uri)
        with self.lock_cache(cache_entry):
            if not cache_entry.exists():
                # open from origin
                with self._open(path_or_uri, "rb") as source, open(
                    cache_entry, "wb"
                ) as cache_source:
                    cache_source.write(source.read())
        return cache_entry

    def _open(self, path_or_uri, mode):
        from smart_open import open

        if str(path_or_uri).startswith("git+file:"):
            return io.BytesIO(git_content(path_or_uri).encode())

        try:
            return open(path_or_uri, mode)
        except Exception as e:
            raise WorkflowError("Failed to open source file {}".format(path_or_uri), e)
