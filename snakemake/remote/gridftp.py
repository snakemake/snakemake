__author__ = "Johannes Köster"
__copyright__ = "Copyright 2017-2019, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"


import os
import subprocess as sp
import time
import shutil

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import os_sync


if not shutil.which("globus-url-copy"):
    raise WorkflowError(
        "The globus-url-copy command has to be available for " "gridftp remote support."
    )

if not shutil.which("gfal-ls"):
    raise WorkflowError(
        "The gfal-* commands need to be available for " "gridftp remote support."
    )


from snakemake.remote import gfal


class RemoteProvider(gfal.RemoteProvider):
    pass


class RemoteObject(gfal.RemoteObject):
    def _globus(self, *args):
        retry = self.provider.retry
        cmd = ["globus-url-copy"] + list(args)
        for i in range(retry + 1):
            try:
                logger.debug(" ".join(cmd))
                return sp.run(
                    cmd, check=True, stderr=sp.PIPE, stdout=sp.PIPE
                ).stdout.decode()
            except sp.CalledProcessError as e:
                if i == retry:
                    raise WorkflowError(
                        "Error calling globus-url-copy:\n{}".format(
                            cmd, e.stderr.decode()
                        )
                    )
                else:
                    # try again after some seconds
                    time.sleep(1)
                    continue

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

            self._globus(
                "-parallel", "4", "-create-dest", "-recurse", "-dp", source, target
            )

            os_sync()
            return self.local_file()
        return None

    def upload(self):
        target = self.remote_file()
        source = "file://" + os.path.abspath(self.local_file())

        if self.exists():
            # first delete file, such that globus does not fail
            self._gfal("rm", target)

        self._globus(
            "-parallel", "4", "-create-dest", "-recurse", "-dp", source, target
        )
