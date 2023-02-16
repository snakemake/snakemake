import time
from snakemake.remote import AbstractRemoteProvider, PooledDomainObject
import snakemake.io
import os
from snakemake.exceptions import (
    WorkflowError
)

class RemoteProvider(AbstractRemoteProvider):
    supports_default = True
    allows_directories = True

    def __init__(
        self, *args, keep_local=False, stay_on_remote=True, is_default=False, **kwargs
    ):
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs
        )

    def remote(self, value, *args, encrypt_data_channel=None, **kwargs):
        return value

    async def inventory(self, cache: snakemake.io.IOCache):
        # for local files, perform BFS via os.scandir to determine existence of files
        if cache.remaining_wait_time <= 0:
            # No more time to create inventory.
            return

        start_time = time.time()

        folders = self.split("/")[:-1]
        if not folders:
            return

        if os.path.isabs(self):
            # For absolute paths, only use scan the immediate parent
            ancestors = [os.path.dirname(self)]
        else:
            ancestors = ["/".join(folders[:i]) for i in range(1, len(folders) + 1)]

        for i, path in enumerate(ancestors):
            if path in cache.exists_local.has_inventory:
                # This path was already scanned before, hence we can stop.
                break
            try:
                with os.scandir(path) as scan:
                    for entry in scan:
                        cache.exists_local[entry.path] = True
                cache.exists_local[path] = True
                cache.exists_local.has_inventory.add(path)
            except FileNotFoundError:
                # Not found, hence, all subfolders cannot be present as well
                for path in ancestors[i:]:
                    cache.exists_local[path] = False
                    cache.exists_local.has_inventory.add(path)
                break
            except PermissionError:
                raise WorkflowError(
                    "Insufficient permissions to access {}. "
                    "Please make sure that all accessed files and directories "
                    "are readable and writable for you.".format(self)
                )

        cache.remaining_wait_time -= time.time() - start_time