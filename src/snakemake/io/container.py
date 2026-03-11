from __future__ import annotations

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import pickle
from typing import (
    Callable,
    List,
)


from snakemake.io import Namedlist


class InputFiles(Namedlist):
    def _predicated_size_files(self, predicate: Callable) -> List[int]:
        async def sizes() -> List[int]:
            async def get_size(f: _IOFile) -> Optional[int]:
                if await predicate(f):
                    return await f.size()
                return None

            sizes = await asyncio.gather(*map(get_size, self))
            return [res for res in sizes if res is not None]

        # the async_run method is expected to be provided through the eval context
        return globals().get("async_run", async_run_fallback)(sizes())

    @property
    def size_files(self):
        async def func_true(_) -> bool:
            return True

        return self._predicated_size_files(func_true)

    @property
    def size_tempfiles(self):
        async def is_temp(iofile):
            return is_flagged(iofile, "temp")

        return self._predicated_size_files(is_temp)

    @property
    def size_files_kb(self):
        return [f / 1024 for f in self.size_files]

    @property
    def size_files_mb(self):
        return [f / 1024 for f in self.size_files_kb]

    @property
    def size_files_gb(self):
        return [f / 1024 for f in self.size_files_mb]

    @property
    def size(self):
        return sum(self.size_files)

    @property
    def size_kb(self):
        return sum(self.size_files_kb)

    @property
    def size_mb(self):
        return sum(self.size_files_mb)

    @property
    def temp_size_mb(self):
        return sum(self.size_tempfiles)

    @property
    def size_gb(self):
        return sum(self.size_files_gb)


class OutputFiles(Namedlist):
    pass


class Wildcards(Namedlist):
    pass


class Params(Namedlist):
    pass


class ResourceList(Namedlist):
    pass


class Log(Namedlist):
    pass
