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


if not shutil.which("gfal-copy"):
    raise WorkflowError("The gfal-* commands need to be available for "
                        "gfal remote support.")


class RemoteProvider(AbstractRemoteProvider):

    supports_default = True
    allows_directories = True

    def __init__(self, *args, stay_on_remote=False, retry=5, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)
        self.retry = retry

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "gsiftp://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        # TODO gfal provides more. Extend this list.
        return ["gsiftp://", "srm://"]


class RemoteObject(AbstractRemoteObject):
    mtime_re = re.compile(r"^\s*Modify: (.+)$", flags=re.MULTILINE)
    size_re = re.compile(r"^\s*Size: ([0-9]+).*$", flags=re.MULTILINE)

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)

    def _gfal(self,
              cmd, *args,
              retry=None,
              raise_workflow_error=True):
        if retry is None:
            retry = self.provider.retry
        _cmd = ["gfal-" + cmd] + list(args)
        for i in range(retry + 1):
            try:
                logger.debug(_cmd)
                return sp.run(_cmd,
                              check=True,
                              stderr=sp.PIPE,
                              stdout=sp.PIPE).stdout.decode()
            except sp.CalledProcessError as e:
                if i == retry:
                    if raise_workflow_error:
                        raise WorkflowError("Error calling gfal-{}:\n{}".format(
                            cmd, e.stderr.decode()))
                    else:
                        raise e
                else:
                    # try again after some seconds
                    time.sleep(1)
                    continue

    # === Implementations of abstract class members ===

    def exists(self):
        try:
            self._gfal("ls", "-a", self.remote_file(),
                       retry=0, raise_workflow_error=False)
        except sp.CalledProcessError as e:
            if e.returncode == 2:
                # exit code 2 means no such file or directory
                return False
            else:
                raise WorkflowError("Error calling gfal-ls:\n{}".format(
                    e.stderr.decode()))
        # exit code 0 means that the file is present
        return True

    def _stat(self):
        stat = self._gfal("stat", self.remote_file())
        return stat

    def mtime(self):
        # assert self.exists()
        stat = self._stat()
        mtime = self.mtime_re.search(stat).group(1)
        date = datetime.strptime(mtime, "%Y-%m-%d %H:%M:%S.%f")
        return date.timestamp()

    def size(self):
        # assert self.exists()
        stat = self._stat()
        size = self.size_re.search(stat).group(1)
        return int(size)

    def download(self):
        if self.exists():
            if self.size() == 0:
                # Globus erroneously thinks that a transfer is incomplete if a
                # file is empty. Hence we manually touch the local file.
                self.local_touch_or_create()
                return self.local_file()
            # Download file. Wait for staging.
            source = self.remote_file()
            target = "file://" + os.path.abspath(self.local_file())

            # disable all timeouts (file transfers can take a long time)
            self._gfal("copy", "-p", "-f", "-n", "4", "-t", "0", "-T", "0",
                       source, target)

            os.sync()
            return self.local_file()
        return None

    def upload(self):
        target = self.remote_file()
        source = "file://" + os.path.abspath(self.local_file())
        # disable all timeouts (file transfers can take a long time)
        self._gfal("copy", "-p", "-f", "-n", "4", "-t", "0", "-T", "0",
                   source, target)

    @property
    def list(self):
        # TODO implement listing of remote files with patterns
        raise NotImplementedError()

    def host(self):
        return self.local_file().split("/")[0]
