__author__ = "Johannes Köster"
__copyright__ = "Copyright 2017, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import os
import re
import shutil
import subprocess as sp
from datetime import datetime

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError
from snakemake.common import lazy_property
from snakemake.logging import logger


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
        cmd = ["uberftp"] + list(args)
        logger.debug(" ".join(cmd))
        try:
            return sp.run(cmd, **kwargs)
        except sp.CalledProcessError as e:
            raise WorkflowError("Error calling uberftp.", e)

    def _uberftp_exists(self, url):
        res = self._uberftp("-ls", url, stdout=sp.PIPE, stderr=sp.STDOUT)
        return "No match for" not in res.stdout.decode()

    # === Implementations of abstract class members ===

    def exists(self):
        return self._uberftp_exists(self.remote_file())

    def mtime(self):
        assert self.exists()
        res = self._uberftp("-ls", self.remote_file(), check=True, stdout=sp.PIPE)
        date = " ".join(res.stdout.decode().split()[-4:-1])
        # first, try minute resolution
        try:
            date = datetime.strptime(date, "%b %d %Y")
        except:
            date = datetime.strptime(date, "%b %d %H:%M")
            date.replace(year=datetime.now().year)
        return date.timestamp()

    def size(self):
        assert self.exists()
        res = self._uberftp("-size", self.remote_file(), check=True, stdout=sp.PIPE)
        return int(res.stdout.decode())

    def download(self):
        if self.exists():
            os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)
            self._uberftp(self.remote_file(),
                          "file://" + os.path.abspath(self.local_file()),
                          check=True)
            os.sync()
            return self.local_file()
        return None

    def upload(self):
        self._uberftp("-rm", self.remote_file())
        prefix = self.protocol[:-1] # remove last slash
        for d in self.local_file().split("/")[:-1]:
            prefix += "/" + d
            if not self._uberftp_exists(prefix):
                self._uberftp("-mkdir", prefix)
        self._uberftp("file://" + os.path.abspath(self.local_file()),
                      self.remote_file(), check=True)

    @property
    def list(self):
        # TODO implement listing of remote files with patterns
        return []
