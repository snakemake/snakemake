__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import hashlib
from pathlib import Path
import re
import os

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
        self.cacheable_prefixes = re.compile("|".join(self.cache_whitelist))

    def lock_cache(self, entry):
        from filelock import FileLock

        return Filelock(entry.with_suffix(".lock"))

    def is_cacheable(self, path_or_uri):
        # TODO remove special git url handling once included in smart_open
        if path_or_uri.startswith("git+file:"):
            return False
        return is_local_file(path_or_uri) and self.cacheable_prefixes.match(path_or_uri)

    def open(self, path_or_uri, mode="r"):
        from smart_open import open

        # TODO add git support to smart_open!
        if self.is_cacheable(path_or_uri):
            # hash the uri and check cache
            urihash = hashlib.sha256()
            urihash.update(path_or_uri.encode())
            urihash = urihash.hexdigest()
            cache_entry = self.cache / urihash

            with self.lock_cache(cache_entry):
                if cache_entry.exists:
                    # open from cache
                    return self._open(cache_entry, mode)
                else:
                    # open from origin
                    with self._open(path_or_uri, "rb") as source, open(
                        cache_entry, "wb"
                    ) as cache_source:
                        cache_source.write(source)
        return self._open(path_or_uri, mode)

    def _open(self, path_or_uri, mode):
        from smart_open import open

        if path_or_uri.startswith("git+file:"):
            return git_content(path_or_uri)

        try:
            return open(path_or_uri, mode)
        except Exception as e:
            raise WorkflowError("Failed to open source file {}".format(path_or_uri), e)
