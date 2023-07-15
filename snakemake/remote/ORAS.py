__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import os
import time

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
import snakemake.io

try:
    import oras
    import oras.provider
except ImportError as e:
    raise WorkflowError(
        "The Python 3 package 'oras` needs to be installed for this remote. %s" % e.msg
    )

# Keep cached oras inventories for a session
caches = {}


class ORASCache(oras.provider.Registry):
    """
    An ORAS Cache is intended to hold one or more data files.

    To properly index, when the artifact is being created, it is recommended
    to upload the files listing individually instead of reference of a
    directory to make it possible to index with Snakemake.
    """

    def exists(self, artifact, filename):
        """
        Determine if a reference exists
        """
        if artifact not in caches:
            self.inventory(artifact, filename)

        if artifact not in caches or filename not in caches[artifact]:
            return False
        return filename in caches[artifact]

    def list(self, filename):
        """
        Given a path
        """

    def download(self, artifact, filename):
        """
        Download a path (registry and filename) to a destination.
        """
        container = self.get_container(artifact)
        digest = caches[artifact][filename]["digest"]
        dest = caches[artifact][filename]["annotations"][
            "org.opencontainers.image.title"
        ]
        self.download_blob(container, digest, dest)

    def size(self, artifact, filename):
        """
        Get the size of path.
        """
        if self.exists(artifact, filename):
            return caches[artifact][filename]["size"]
        return 0

    def parent(self, filename):
        """
        Derive the parent of a path.
        """
        return os.path.dirname(filename)

    def inventory(self, artifact, filename):
        """
        Get the inventory of a registry.
        """
        container = self.get_container(artifact)

        # Get the manifest with file contents
        manifest = self.get_manifest(container)

        # Organize by path. ORAS stories the path (name of file) under org.opencontainers.image.title
        contents = {}
        for layer in manifest.get("layers", {}):
            if "org.opencontainers.image.title" not in layer.get("annotations", {}):
                logger.warning(
                    f"Layer {layer} is missing ORAS annotation org.opencontainers.image.title and cannot be used."
                )
                continue
            path = layer["annotations"]["org.opencontainers.image.title"]
            contents[path] = layer
        caches[artifact] = contents
        return contents


class RemoteProvider(AbstractRemoteProvider):
    supports_default = True

    def __init__(
        self, *args, keep_local=False, stay_on_remote=False, is_default=False, **kwargs
    ):
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs,
        )
        self.client = ORASCache()

    def remote_interface(self):
        return self.client

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "oras://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["oras://"]


class RemoteObject(AbstractRemoteObject):
    def __init__(
        self, *args, keep_local=False, provider=None, user_project=None, **kwargs
    ):
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )
        if "artifact" not in kwargs:
            raise WorkflowError(
                "The ORAS provider requires an artifact URI from an OCI registry."
            )

        if provider:
            self.client = provider.remote_interface()
        else:
            self.client = ORASCache()

        # keep user_project available for when bucket is initialized
        self._user_project = user_project
        self._artifact = kwargs["artifact"]

    async def inventory(self, cache: snakemake.io.IOCache):
        if cache.remaining_wait_time <= 0:
            # No more time to create inventory.
            return

        contents = self.client.inventory(self._artifact, self._file)

        start_time = time.time()
        for name, metadata in contents.items():
            cache.exists_remote[name] = True
            cache.size[name] = metadata["size"]
        cache.remaining_wait_time -= time.time() - start_time

    def get_inventory_parent(self):
        return self.client.parent(self._file)

    def exists(self):
        return self.client.exists(self._artifact, self._file)

    def mtime(self):
        if self.client.exists(self._artifact, self._file):
            return 0
        else:
            raise WorkflowError(
                "The file does not seem to exist remotely: %s" % self.local_file()
            )

    def size(self):
        size = self.client.size(self._artifact, self._file)
        if size:
            return size

    def _download(self):
        """Download with maximum retry duration of 600 seconds (10 minutes)"""
        if not self.exists():
            return None
        return self.client.download(self._artifact, self._file)

    def is_directory(self):
        # TODO not sure if an empty path can be represented in an archive
        if snakemake.io.is_flagged(self.file(), "directory"):
            return True
        elif self.client.exists(self._artifact, self.file()):
            return False
        else:
            return any(self.directory_entries())

    @property
    def name(self):
        return self._file

    @property
    def list(self):
        return self.client.list(self._file)
