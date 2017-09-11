__author__ = "Oliver Stolpe"
__copyright__ = "Copyright 2017, BIH Core Unit Bioinformatics"
__email__ = "oliver.stolpe@bihealth.org"
__license__ = "MIT"

import os
from contextlib import contextmanager
from datetime import datetime, timedelta
from pytz import timezone

# module-specific
from snakemake.remote import AbstractRemoteProvider, DomainObject
from snakemake.exceptions import iRODSFileException, WorkflowError

try:
    # third-party modules
    from irods.session import iRODSSession
    from irods.meta import iRODSMeta
    from irods.models import DataObject
    import irods.keywords as kw
except ImportError as e:
    raise WorkflowError("The Python 3 package 'python.irods' " +
        "must be installed to use iRODS remote() file functionality. %s" % e.msg)


utc = datetime.utcfromtimestamp(0)


class RemoteProvider(AbstractRemoteProvider):

    supports_default = True

    def __init__(self, *args, stay_on_remote=False, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return ''

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return []


class RemoteObject(DomainObject):
    """ This is a class to interact with an iRODS server.
    """

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)

    # === Implementations of abstract class members ===

    @contextmanager #makes this a context manager. after 'yield' is __exit__()
    def irods_session(self):
        # if args have been provided to remote(), use them over those given to RemoteProvider()
        args_to_use = self.provider.args
        if len(self.args):
            args_to_use = self.args

        zone = self.remote_path.split(os.sep)[1]

        # use kwargs passed in to remote() to override those given to the RemoteProvider()
        # default to the host and port given as part of the file, falling back to one specified
        # as a kwarg to remote() or the RemoteProvider (overriding the latter with the former if both)
        kwargs_to_use = {}
        kwargs_to_use["host"] = self.host
        kwargs_to_use["port"] = int(self.port) if self.port else 1247
        kwargs_to_use["zone"] = self.provider.kwargs.get('zone', zone)
        kwargs_to_use["user"] = self.provider.kwargs['user']
        kwargs_to_use["password"] = self.provider.kwargs['password']

        # do no take everything from provider.kwargs because then it also passes the timezone
        for k,v in self.kwargs.items():
            # allow only for the previously defined parameters
            # otherwise the 'overwrite' parameter would also be passed
            if k not in kwargs_to_use.keys():
                continue
            kwargs_to_use[k] = v

        session = iRODSSession(**kwargs_to_use)
        yield session
        session.cleanup()

    def exists(self):
        if self._matched_address:
            with self.irods_session() as session:
                try:
                    obj = session.data_objects.get(self.remote_path)
                    return True
                except:
                    return False
        else:
            raise iRODSFileException("The file cannot be parsed as an iRODS path in form 'host:port/path/to/file': %s" % self.local_file())

    def _convert_time(self, timestamp, tz=None):
        dt = timestamp.replace(tzinfo=timezone('UTC'))
        if tz:
            dt = dt.astimezone(timezone(tz))
        return dt

    def mtime(self):
        if self.exists():
            with self.irods_session() as session:
                obj = session.data_objects.get(self.remote_path)
                meta = session.metadata.get(DataObject, self.remote_path)

                # if mtime was set in metadata during upload, take this information
                # otherwise fall back to iRODS timestamp (upload time!) and change
                # timezone accordingly (iRODS might ignore the servers local timezone)
                for m in meta:
                    if m.name == 'mtime':
                        mtime = float(m.value)
                        break
                else:
                    dt = self._convert_time(obj.modify_time, self.provider.kwargs.get('timezone', None))
                    utc2 = self._convert_time(utc)
                    mtime = (dt - utc2).total_seconds()
                
                return int(mtime)
        else:
            raise iRODSFileException("The file does not seem to exist remotely: %s" % self.local_file())

    def atime(self):
        if self.exists():
            with self.irods_session() as session:
                obj = session.data_objects.get(self.remote_path)
                meta = session.metadata.get(DataObject, self.remote_path)

                # if mtime was set in metadata during upload, take this information
                # otherwise fall back to iRODS timestamp (upload time!) and change
                # timezone accordingly (iRODS might ignore the servers local timezone)
                for m in meta:
                    if m.name == 'atime':
                        atime = float(m.value)
                        break
                else:
                    dt = self._convert_time(obj.modify_time, self.provider.kwargs.get('timezone', None))
                    utc2 = self._convert_time(utc)
                    atime = (dt - utc2).total_seconds()
                
                return int(atime)
        else:
            raise iRODSFileException("File doesn't exist remotely: %s" % self.local_file())

    def is_newer(self, time):
        """ Returns true of the file is newer than time, or if it is
            a symlink that points to a file newer than time. """
        return self.mtime() > time
        
    def size(self):
        if self.exists():
            with self.irods_session() as session:
                obj = session.data_objects.get(self.remote_path)
                return int(obj.size)
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        with self.irods_session() as session:
            if self.exists():
                if make_dest_dirs:
                    os.makedirs(os.path.dirname(self.local_path), exist_ok=True)

                # force irods to overwrite existing file if this option is set
                opt = {}
                if self.kwargs.get('overwrite', None):
                    opt[kw.FORCE_FLAG_KW] = ''

                # get object and change timestamp
                obj = session.data_objects.get(self.remote_path, self.local_path, options=opt)
                os.utime(self.local_path, (self.atime(), self.mtime()))
                os.sync()
            else:
                raise iRODSFileException("The file does not seem to exist remotely: %s" % self.local_file())

    def upload(self):
        with self.irods_session() as session:
            # get current local timestamp
            stat = os.stat(self.local_path)

            # create folder structure on remote
            folders = os.path.dirname(self.remote_path).split(os.sep)[1:]
            collpath = os.sep

            for folder in folders:
                collpath = os.path.join(collpath, folder)

                try:
                    print("trying to get {}".format(collpath))
                    session.collections.get(collpath)
                except:
                    print("creating {}".format(collpath))
                    session.collections.create(collpath)

            # upload file and store local timestamp in metadata since irods sets the files modification time to
            # the upload time rather than retaining the local modification time
            session.data_objects.put(self.local_path, self.remote_path)
            session.metadata.add(DataObject, self.remote_path, iRODSMeta('mtime', str(stat.st_mtime), 's'))
            session.metadata.add(DataObject, self.remote_path, iRODSMeta('atime', str(stat.st_atime), 's'))
            session.metadata.add(DataObject, self.remote_path, iRODSMeta('ctime', str(stat.st_ctime), 's'))

    @property
    def list(self):
        file_list = []

        first_wildcard = self._iofile.constant_prefix()
        dirname = os.path.dirname(first_wildcard.replace(self.path_prefix, ""))

        with self.irods_session() as session:
            collection = session.collections.get(dirname)
            for current_collection, subcollections, objs in collection.walk():
                for obj in objs:
                    file_list.append(obj.path.lstrip('/'))

        return file_list
