__author__ = "Johannes Köster"
__copyright__ = "Copyright 2017, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"


import subprocess as sp

from snakemake.remote import gfal


class RemoteProvider(gfal.RemoteProvider):
    pass


class RemoteObject(gfal.RemoteObject):
    def _globus(self,
                *args):
        retry = self.provider.retry
        cmd = ["globus-url-copy"] + list(args)
        for i in range(retry + 1):
            try:
                logger.debug(" ".join(cmd))
                return sp.run(cmd,
                              check=True,
                              stderr=sp.PIPE,
                              stdout=sp.PIPE).stdout.decode()
            except sp.CalledProcessError as e:
                if i == retry:
                    raise WorkflowError("Error calling gfal-{}:\n{}".format(
                        cmd, e.stderr.decode()))
                else:
                    # try again after some seconds
                    time.sleep(1)
                    continue

    def download(self):
        if self.exists():
            # Download file. Wait for staging.
            source = self.remote_file()
            target = "file://" + os.path.abspath(self.local_file())

            self._globus("-parallel", "4", "-create-dest", "-recurse",
                         source, target)

            os.sync()
            return self.local_file()
        return None

    def upload(self):
        target = self.remote_file()
        source = "file://" + os.path.abspath(self.local_file())

        self._globus("-parallel", "4", "-create-dest", "-recurse",
                     source, target)
