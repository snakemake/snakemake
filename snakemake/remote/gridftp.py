__author__ = "Johannes Köster"
__copyright__ = "Copyright 2017, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import os
import re
import shutil
import subprocess as sp

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError
from snakemake.common import lazy_property


if not shutil.which("uberftp"):
    raise WorkflowError("The uberftp command needs to be available for GridFTP support")


class RemoteProvider(AbstractRemoteProvider):

    supports_default = True

    def __init__(self, *args, stay_on_remote=False, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "gsiftp://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["gsiftp://"]


class RemoteObject(AbstractRemoteObject):

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)

    def _uberftp(self, *args, **kwargs):
        try:
            sp.run(["uberftp"] + args, check=True)
        except sp.CalledProcessError as e:
            raise WorkflowError("Error calling uberftp.", e)

    # === Implementations of abstract class members ===

    def exists(self):
        res = self._uberftp("-ls", self.remote_file(), stdout=sp.PIPE)
        return "No match for" not in res.stdout.decode()

    def mtime(self):
        assert self.exists()
        res = self._uberftp("-ls", self.remote_file(), stdout=sp.PIPE)
        date = " ".join(res.stdout.decode().split()[-4:-1])
        # first, try minute resolution
        try:
            date = datetime.strptime(date, "%b %d %H:%M")
        except:
            date = datetime.strptime(date, "%b %d %Y")

        return date.timestamp()

    def size(self):
        assert self.exists()
        res = self._uberftp("-size", self.remote_file(), stdout=sp.PIPE)
        return int(res.stdout.decode())

    def download(self):
        if self.exists():
            os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)
            self._uberftp(self.remote_file(), "file://" + self.local_file())
            os.sync()
            return self.local_file()
        return None

    def upload(self):
        self._uberftp("-rm", self.remote_file())
        self._uberftp("file://" + self.local_file(), self.remote_file())

    @property
    def list(self):
        # TODO implement listing of remote files with patterns
        return []
