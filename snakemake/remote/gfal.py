__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import os
import stat

from snakemake.remote import (
    AbstractRemoteProvider,
    AbstractRemoteRetryObject,
    check_deprecated_retry,
)
from snakemake.exceptions import WorkflowError
from snakemake.utils import os_sync

try:
    # third-party modules
    import gfal2
except ImportError as e:
    raise WorkflowError(
        'gfal2 python binding not found. The conda package name is "python-gfal2" '
    ) from e


class RemoteProvider(AbstractRemoteProvider):
    supports_default = True

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
    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )
        if self.provider:
            self.gfalcntxt = provider.remote_interface()
        else:
            self.gfalcntxt = gfal2.creat_context()

    # === Implementations of abstract class members ===

    def exists(self):
        try:
            self.gfalcntxt.stat(self.remote_file())
        except gfal2.GError as e:
            if "No such file or directory" in e.message:
                return False
            else:
                raise WorkflowError(
                    f"Error calling gfal2-python stat:\n\t{e.message} \n\t{e.code}"
                ) from e

        return True

    def _stat(self):
        return self._gfal("stat", self.remote_file())

    def mtime(self):
        return self.gfalcntxt.stat(self.remote_file()).st_mtime

    def size(self):
        stat = self.gfalcntxt.stat(self.remote_file())
        return stat.st_size

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
        target = f"file://{os.path.abspath(self.local_file())}"

        self._transfer_file(target, source)

        os_sync()
        return self.local_file()

    def _upload(self):
        target = self.remote_file()
        source = f"file://{os.path.abspath(self.local_file())}"

        self._transfer_file(target, source)

    def _transfer_file(self, target, source):
        params = self.gfalcntxt.transfer_parameters()
        params.overwrite = True
        params.nbstreams = 4
        params.timeout = 1200
        params.create_parent = True

        try:
            self.gfalcntxt.filecopy(params, source, target)
        except gfal2.GError as e:
            raise WorkflowError(
                f"Error calling gfal2-python copy:\n\t{e.message} \n\t{e.code}"
            ) from e

    def remove(self):
        self._do_rm(self.remote_file())

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
        [self._do_rm(base_dir + c) for c in contents if c not in (".", "..")]
        self.gfalcntx.rmdir(surl)

    @property
    def list(self):
        # TODO implement listing of remote files with patterns
        raise NotImplementedError()

    def host(self):
        return self.local_file().split("/")[0]
