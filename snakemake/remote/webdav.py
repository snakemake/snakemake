__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2017, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os, sys
import email.utils
from contextlib import contextmanager
import functools

# module-specific
from snakemake.remote import AbstractRemoteProvider, AbstractRemoteObject, DomainObject
from snakemake.exceptions import WebDAVFileException, WorkflowError
from snakemake.utils import os_sync

try:
    # third-party modules
    import aioeasywebdav
    import asyncio
except ImportError as e:
    raise WorkflowError(
        "The Python 3 packages 'aioeasywebdav' "
        " and 'asyncio' must be present to use WebDAV remote() file "
        "functionality. %s" % e.msg
    )


class RemoteProvider(AbstractRemoteProvider):
    def __init__(
        self, *args, keep_local=False, stay_on_remote=False, is_default=False, **kwargs
    ):
        # loop = asyncio.get_event_loop()
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs
        )

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "https://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["http://", "https://"]


class RemoteObject(DomainObject):
    """ This is a class to interact with a WebDAV file store.
    """

    def __init__(self, *args, keep_local=False, **kwargs):
        # self.loop = asyncio.get_event_loop()
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, **kwargs)

    @contextmanager
    def webdavc(self):
        newloop = False
        if not hasattr(self, "loop"):
            try:
                self.loop = asyncio.get_event_loop()
                if self.loop.is_running():
                    raise NotImplementedError(
                        "Cannot use aioutils in " "asynchroneous environment"
                    )
            except:
                newloop = True
                self.loop = asyncio.new_event_loop()
                asyncio.set_event_loop(self.loop)

            self.loop = asyncio.new_event_loop()
            asyncio.set_event_loop(self.loop)

        if not hasattr(self, "conn") or (
            hasattr(self, "conn") and not isinstance(self.conn, aioeasywebdav.Client)
        ):
            # if args have been provided to remote(), use them over those given to RemoteProvider()
            args_to_use = self.provider.args
            if len(self.args):
                args_to_use = self.args

            # use kwargs passed in to remote() to override those given to the RemoteProvider()
            # default to the host and port given as part of the file, falling back to one specified
            # as a kwarg to remote() or the RemoteProvider (overriding the latter with the former if both)
            kwargs_to_use = {}
            kwargs_to_use["host"] = self.host
            kwargs_to_use["protocol"] = self.protocol
            kwargs_to_use["port"] = int(self.port) if self.port != None else 443
            for k, v in self.provider.kwargs.items():
                kwargs_to_use[k] = v
            for k, v in self.kwargs.items():
                kwargs_to_use[k] = v

            # easywebdav wants the protocol without "://"
            kwargs_to_use["protocol"] = kwargs_to_use["protocol"].replace("://", "")

            # monkey patch aioeasywebdav to noop _rate_calc()
            # since we don't care about download progress and
            # the parent (connection) object may be removed before the
            # sleep coroutine has a chance to be scheduled/finish,
            # and aioeasywebdav only calls close() on __del__()
            async def noop(_):
                pass

            aioeasywebdav.Client._rate_calc = noop

            self.conn = aioeasywebdav.connect(*args_to_use, **kwargs_to_use)
        yield

    # === Implementations of abstract class members ===

    def exists(self):
        with self.webdavc() as webdavc:
            path_to_try = self.webdav_file
            return self.loop.run_until_complete(self.conn.exists(self.webdav_file))

    def mtime(self):
        if self.exists():
            with self.webdavc() as webdavc:
                metadata = self.loop.run_until_complete(
                    self.conn.ls(remote_path=self.webdav_file)
                )[0]
                parsed_date = email.utils.parsedate_tz(metadata.mtime)
                epoch_time = email.utils.mktime_tz(parsed_date)
                return epoch_time
        else:
            raise WorkflowError(
                "The file does not seem to exist remotely: %s" % self.webdav_file
            )

    def size(self):
        if self.exists():
            with self.webdavc() as webdavc:
                metadata = self.loop.run_until_complete(
                    self.conn.ls(remote_path=self.webdav_file)
                )[0]
                return int(metadata.size)
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        if self.exists():
            # if the destination path does not exist, make it
            if make_dest_dirs:
                os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)
            with self.webdavc() as webdavc:
                self.loop.run_until_complete(
                    self.conn.download(self.webdav_file, self.local_file())
                )
                os_sync()  # ensure flush to disk
        else:
            raise WorkflowError(
                "The file does not seem to exist remotely: %s" % self.webdav_file
            )

    def upload(self):
        # make containing folder
        with self.webdavc() as webdavc:
            self.loop.run_until_complete(
                self.conn.mkdirs(os.path.dirname(self.webdav_file))
            )
            self.loop.run_until_complete(
                self.conn.upload(self.local_file(), self.webdav_file)
            )

    @property
    def webdav_file(self):
        filepath = (
            self.local_file().replace(self.host, "").replace(":" + str(self.port), "")
        )
        filepath = filepath if not filepath.startswith("/") else filepath[1:]
        return filepath

    @property
    def name(self):
        return self.local_file()

    @property
    def list(self):
        file_list = []

        first_wildcard = (
            self._iofile.constant_prefix()
            .replace(self.host, "")
            .replace(":" + str(self.port), "")
        )
        dirname = (
            first_wildcard[1:] if first_wildcard.startswith("/") else first_wildcard
        )

        while "//" in dirname:
            dirname = dirname.replace("//", "/")
        dirname = dirname.rstrip("/") + "/"

        with self.webdavc() as webdavc:
            for item in self.loop.run_until_complete(self.conn.ls(dirname)):
                file_list.append(
                    os.path.join(os.path.dirname(dirname), item.name.lstrip("/"))
                )
                file_list.append(
                    os.path.join(
                        self._iofile.constant_prefix(), os.path.basename(item.name)
                    )
                )

        return file_list
