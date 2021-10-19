__author__ = "Oliver Stolpe"
__copyright__ = "Copyright 2017, BIH Core Unit Bioinformatics"
__email__ = "oliver.stolpe@bihealth.org"
__license__ = "MIT"

import os
import re

from contextlib import contextmanager
from datetime import datetime, timedelta
from pytz import timezone

# module-specific
from snakemake.remote import AbstractRemoteProvider, AbstractRemoteObject
from snakemake.exceptions import WorkflowError
from snakemake.utils import os_sync

try:
    # third-party modules
    from irods.session import iRODSSession
    from irods.meta import iRODSMeta
    from irods.models import DataObject
    from irods.exception import CollectionDoesNotExist, DataObjectDoesNotExist
    import irods.keywords as kw
except ImportError as e:
    raise WorkflowError(
        "The Python 3 package 'python-irodsclient' "
        + "must be installed to use iRODS remote() file functionality. %s" % e.msg
    )


utc = datetime.utcfromtimestamp(0)


def _irods_session(*args, **kwargs):
    try:
        irods_env_file = os.environ["IRODS_ENVIRONMENT_FILE"]
    except KeyError:
        irods_env_file = kwargs.get(
            "irods_env_file", os.path.expanduser("~/.irods/irods_environment.json")
        )

    if not os.path.isfile(irods_env_file):
        raise WorkflowError(
            "Expecting iRODS configuration file in %s, but none found."
            % kwargs["irods_env_file"]
        )

    with iRODSSession(irods_env_file=irods_env_file) as session:
        try:
            session.collections.get(
                os.path.join(os.sep, session.zone, "home", session.username)
            )
        except ConnectionError as e:
            raise e

        return session


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
            **kwargs
        )
        self._irods_session = _irods_session(*args, **kwargs)

    def remote_interface(self):
        return self._irods_session

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "irods://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["irods://"]


class RemoteObject(AbstractRemoteObject):
    """This is a class to interact with an iRODS server."""

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )
        if provider:
            self._irods_session = provider.remote_interface()
            self._timezone = provider.kwargs.get("timezone")
        else:
            self._irods_session = _irods_session(*args, **kwargs)

    def exists(self):
        try:
            self._irods_session.data_objects.get(self.remote_path)
            return True
        except (CollectionDoesNotExist, DataObjectDoesNotExist):
            return False

    def _convert_time(self, timestamp, tz=None):
        dt = timestamp.replace(tzinfo=timezone("UTC"))
        if tz:
            dt = dt.astimezone(timezone(tz))
        return dt

    def mtime(self):
        if self.exists():
            obj = self._irods_session.data_objects.get(self.remote_path)
            meta = self._irods_session.metadata.get(DataObject, self.remote_path)

            # if mtime was set in metadata during upload, take this information
            # otherwise fall back to iRODS timestamp (upload time!) and change
            # timezone accordingly (iRODS might ignore the servers local timezone)
            for m in meta:
                if m.name == "mtime":
                    mtime = float(m.value)
                    break
            else:
                dt = self._convert_time(obj.modify_time, self._timezone)
                utc2 = self._convert_time(utc)
                mtime = (dt - utc2).total_seconds()

            return int(mtime)
        else:
            raise WorkflowError(
                "The file does not seem to exist remotely: %s" % self.local_file()
            )

    def atime(self):
        if self.exists():
            obj = self._irods_session.data_objects.get(self.remote_path)
            meta = self._irods_session.metadata.get(DataObject, self.remote_path)

            # if mtime was set in metadata during upload, take this information
            # otherwise fall back to iRODS timestamp (upload time!) and change
            # timezone accordingly (iRODS might ignore the servers local timezone)
            for m in meta:
                if m.name == "atime":
                    atime = float(m.value)
                    break
            else:
                dt = self._convert_time(obj.modify_time, self._timezone)
                utc2 = self._convert_time(utc)
                atime = (dt - utc2).total_seconds()

            return int(atime)
        else:
            raise WorkflowError("File doesn't exist remotely: %s" % self.local_file())

    def is_newer(self, time):
        """Returns true of the file is newer than time, or if it is
        a symlink that points to a file newer than time."""
        return self.mtime() > time

    def size(self):
        if self.exists():
            obj = self._irods_session.data_objects.get(self.remote_path)
            return int(obj.size)
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        if self.exists():
            if make_dest_dirs:
                os.makedirs(os.path.dirname(self.local_path), exist_ok=True)

            # force irods to overwrite existing file if this option is set
            opt = {}
            if self.kwargs.get("overwrite"):
                opt[kw.FORCE_FLAG_KW] = ""

            # get object and change timestamp
            obj = self._irods_session.data_objects.get(
                self.remote_path, self.local_path, options=opt
            )
            os.utime(self.local_path, (self.atime(), self.mtime()))
            os_sync()
        else:
            raise WorkflowError(
                "The file does not seem to exist remotely: %s" % self.local_file()
            )

    def upload(self):
        # get current local timestamp
        stat = os.stat(self.local_path)

        # create folder structure on remote
        folders = os.path.dirname(self.remote_path).split(os.sep)[1:]
        collpath = os.sep

        for folder in folders:
            collpath = os.path.join(collpath, folder)

            try:
                self._irods_session.collections.get(collpath)
            except:
                self._irods_session.collections.create(collpath)

        # upload file and store local timestamp in metadata since irods sets the files modification time to
        # the upload time rather than retaining the local modification time
        self._irods_session.data_objects.put(self.local_path, self.remote_path)

        # erase meta data (if exists) before adding it. there is no update routine available in the API
        for m in self._irods_session.metadata.get(DataObject, self.remote_path):
            if m.name in ("mtime", "atime", "ctime"):
                self._irods_session.metadata.remove(
                    DataObject, self.remote_path, iRODSMeta(m.name, m.value, m.units)
                )

        self._irods_session.metadata.add(
            DataObject, self.remote_path, iRODSMeta("mtime", str(stat.st_mtime), "s")
        )
        self._irods_session.metadata.add(
            DataObject, self.remote_path, iRODSMeta("atime", str(stat.st_atime), "s")
        )
        self._irods_session.metadata.add(
            DataObject, self.remote_path, iRODSMeta("ctime", str(stat.st_ctime), "s")
        )

    @property
    def remote_path(self):
        return os.path.join(os.sep, self._irods_session.zone, self.local_file())

    @property
    def name(self):
        return self.local_file()

    @property
    def local_path(self):
        return self.local_file()

    @property
    def list(self):
        file_list = []

        first_wildcard = self._iofile.constant_prefix()
        dirname = os.path.dirname(first_wildcard)

        collection = self._irods_session.collections.get(dirname)
        for _, _, objs in collection.walk():
            for obj in objs:
                file_list.append(obj.path.lstrip("/"))

        return file_list
