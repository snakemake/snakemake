__author__ = "Johannes Köster"
__copyright__ = "Copyright 2017, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import os
import re
import shutil
import subprocess as sp
from datetime import datetime
import time

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError
from snakemake.common import lazy_property
from snakemake.logging import logger


if not shutil.which("uberftp"):
    raise WorkflowError("The uberftp command needs to be available for "
                        "GridFTP support.")
if not shutil.which("globus-url-copy"):
    raise WorkflowError("The globus-url-copy command needs to be available for "
                        "GridFTP support.")



class RemoteProvider(AbstractRemoteProvider):

    supports_default = True

    def __init__(self, *args, stay_on_remote=False, retry=10, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)
        self.retry = retry

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

    def _globus_url_copy(self, source, target):
        try:
            return sp.run(["globus-url-copy", "-fast", "-cd", "-r",
                          "-rst", "-rst-retries", self.provider.retry,
                          source, target], check=True, stderr=sp.PIPE)
        except sp.CalledProcessError as e:
            raise WorkflowError("Error calling globus-url-copy:\n{}".format(
                e.stderr.decode()))

    def _uberftp(self, *args, **kwargs):
        if "stderr" not in kwargs:
            kwargs["stderr"] = sp.PIPE
        # try several times for error robustness
        for i in range(self.provider.retry + 1):
            cmd = " ".join(["uberftp"] + list(args))
            logger.debug(cmd)
            try:
                # use shell=True because otherwise login seeems to be
                # unreliable
                return sp.run(cmd, **kwargs, shell=True)
            except sp.CalledProcessError as e:
                if i == 9:
                    raise WorkflowError("Error calling uberftp:\n{}".format(
                        e.stderr.decode()))
                time.sleep(0.5)

    def _uberftp_exists(self, url):
        for i in range(self.provider.retry + 1):
            res = self._uberftp("-ls", url, stdout=sp.PIPE, stderr=sp.STDOUT)
            if res.returncode == 0:
                return True
            elif "No match for" in res.stdout.decode():
                return False
            else:
                if i == 9:
                    raise WorkflowError(
                        "Error calling uberftp:\n{}".format(
                            res.stdout.decode()))
                time.sleep(0.5)

    # === Implementations of abstract class members ===

    def exists(self):
        return self._uberftp_exists(self.remote_file())

    def mtime(self):
        # assert self.exists()
        res = self._uberftp("-ls", self.remote_file(), check=True, stdout=sp.PIPE)
        date = " ".join(res.stdout.decode().split()[-4:-1])
        # first, try minute resolution
        try:
            date = datetime.strptime(date, "%b %d %Y")
        except:
            date = datetime.strptime(date, "%b %d %H:%M")
            date = date.replace(year=datetime.now().year)
        return date.timestamp()

    def size(self):
        # assert self.exists()
        res = self._uberftp("-size", self.remote_file(), check=True, stdout=sp.PIPE)
        return int(res.stdout.decode())

    def download(self):
        if self.exists():
            os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)

            # Download file. Wait for staging.
            source = self.remote_file()
            target = "file://" + os.path.abspath(self.local_file())

            self._globus_url_copy(source, target)

            os.sync()
            return self.local_file()
        return None

    def upload(self):
        if os.path.isdir(self.local_file()):
            raise WorkflowError("Directories are not supported by gridftp remote.")
        if self.exists():
            self._uberftp("-rm", self.remote_file(), check=True)

        prefix = self.protocol + self.host()
        # omit first and last elements (host and file)
        for d in self.local_file().split("/")[1:-1]:
            prefix += "/" + d
            if not self._uberftp_exists(prefix):
                self._uberftp("-mkdir", prefix, check=True)

        # Upload file.
        source = "file://" + os.path.abspath(self.local_file())
        target = self.remote_file()

        self._globus_url_copy(source, target)

    @property
    def list(self):
        # TODO implement listing of remote files with patterns
        return []

    def host(self):
        return self.local_file().split("/")[0]
