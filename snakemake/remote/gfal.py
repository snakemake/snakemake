__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import os
import re
import stat
import shutil
import subprocess as sp
from datetime import datetime

from snakemake.remote import (
    AbstractRemoteObject,
    AbstractRemoteProvider,
    AbstractRemoteRetryObject,
    check_deprecated_retry,
)
from snakemake.exceptions import WorkflowError
from snakemake.common import lazy_property
from snakemake.logging import logger
from snakemake.utils import os_sync

gfal_python = True
try:
    # third-party modules
    import gfal2
except ImportError as e:
    gfal_python = False
    logger.warning(
        "gfal2 python binding not found. Will use slower gfal2-* commands as fallback."
    )
    if not shutil.which("gfal-copy"):
        raise WorkflowError(
            "Gfal-* are not found. Gfal-* commands are needed to be available for gfal remote support."
        )


class RemoteProvider(AbstractRemoteProvider):
    supports_default = True
    allows_directories = True
    if gfal_python:
        gfalcntx = gfal2.creat_context()

    def __init__(
        self,
        *args,
        keep_local=False,
        stay_on_remote=False,
        is_default=False,
        retry=None,
        **kwargs,
    ):
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs,
        )
        check_deprecated_retry(retry)

    def remote_interface(self):
        return self.gfalcntx

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "gsiftp://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        # TODO gfal provides more. Extend this list.
        return ["gsiftp://", "srm://"]


class RemoteObject(AbstractRemoteRetryObject):

    mtime_re = re.compile(r"^\s*Modify: (.+)$", flags=re.MULTILINE)
    size_re = re.compile(r"^\s*Size: ([0-9]+).*$", flags=re.MULTILINE)

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )
        if gfal_python:
            if self.provider:
                self.gfalcntxt = provider.remote_interface()
            else:
                self.gfalcntxt = gfal2.creat_context()

    def _gfal(self, cmd, *args, retry=None, raise_workflow_error=True):
        check_deprecated_retry(retry)
        _cmd = ["gfal-" + cmd] + list(args)
        try:
            logger.debug(_cmd)
            return sp.run(
                _cmd, check=True, stderr=sp.PIPE, stdout=sp.PIPE
            ).stdout.decode()
        except sp.CalledProcessError as e:
            if raise_workflow_error:
                raise WorkflowError(
                    f"Error calling gfal-{cmd}:\n{e.stderr.decode()}"
                ) from e
            else:
                raise e

    # === Implementations of abstract class members ===

    def exists(self):

        if gfal_python:
            try:
                self.gfalcntxt.stat(self.remote_file())
            except gfal2.GError as e:
                if "No such file or directory" in e.message:
                    return False
                else:
                    raise WorkflowError(
                        f"Error calling gfal2-python stat:\n\t{e.message} \n\t{e.code}"
                    ) from e

        else:
            try:
                self._gfal("ls", "-a", self.remote_file(), raise_workflow_error=False)
            except sp.CalledProcessError as e:
                if e.returncode == 2:
                    # exit code 2 means no such file or directory
                    return False
                else:
                    raise WorkflowError(
                        f"Error calling gfal-ls:\n{e.stderr.decode()}"
                    ) from e

        return True

    def _stat(self):
        return self._gfal("stat", self.remote_file())

    def mtime(self):
        if gfal_python:
            return self.gfalcntxt.stat(self.remote_file()).st_mtime

        # assert self.exists()
        stat = self._stat()
        mtime = self.mtime_re.search(stat).group(1)
        date = datetime.strptime(mtime, "%Y-%m-%d %H:%M:%S.%f")
        return date.timestamp()

    def size(self):
        if gfal_python:
            stat = self.gfalcntxt.stat(self.remote_file())
            return stat.st_size

        # assert self.exists()
        stat = self._stat()
        size = self.size_re.search(stat).group(1)
        return int(size)

    def _download(self):
        if not self.exists():
            return None

        if self.size() == 0:
            # Globus erroneously thinks that a transfer is incomplete if a
            # file is empty. Hence we manually touch the local file.
            self.local_touch_or_create()
            return self.local_file()
        # Download file. Wait for staging.
        source = self.remote_file()
        target = "file://" + os.path.abspath(self.local_file())

        if gfal_python:
            self._transfer_file(target, source)
        else:
            # disable all timeouts (file transfers can take a long time)
            self._gfal(
                "copy", "-p", "-f", "-n", "4", "-t", "0", "-T", "0", source, target
            )

        os_sync()
        return self.local_file()

    def _upload(self):
        target = self.remote_file()
        source = "file://" + os.path.abspath(self.local_file())
        if gfal_python:
            self._transfer_file(target, source)
        else:
            # disable all timeouts (file transfers can take a long time)
            self._gfal(
                "copy", "-p", "-f", "-n", "4", "-t", "0", "-T", "0", source, target
            )

    def _transfer_file(self, target, source):
        params = self.gfalcntxt.transfer_parameters()
        params.overwrite = True
        params.nbstreams = 4
        params.timeout = 1200
        params.create_parent = True

        try:
            r = self.gfalcntxt.filecopy(params, source, target)
        except gfal2.GError as e:
            raise WorkflowError(
                f"Error calling gfal2-python copy:\n\t{e.message} \n\t{e.code}"
            ) from e

    def remove(self):
        if gfal_python:
            self._do_rm(self.remote_file())
        self._gfal("rm", "-r", self.remote_file())

    def _do_rm(self, surl):
        """ """
        st = self.gfalcntxt.stat(surl)
        if stat.S_ISDIR(st.st_mode):
            self._do_rmdir(surl)
            return
        self.gfalcntxt.unlink(surl)

    def _do_rmdir(self, surl):
        """ """
        base_dir = surl
        if base_dir.endswith("/"):
            base_dir = f"{base_dir}/"
        contents = self.gfalcntxt.listdir(surl)
        [self._do_rm(base_dir + c) for c in contents if c not in [".", ".."]]
        self.gfalcntx.rmdir(surl)

    @property
    def list(self):
        # TODO implement listing of remote files with patterns
        raise NotImplementedError()

    def host(self):
        return self.local_file().split("/")[0]
