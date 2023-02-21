__author__ = "Chris Burr"
__copyright__ = "Copyright 2022, Chris Burr"
__email__ = "christopher.burr@cern.ch"
__license__ = "MIT"

import os
from os.path import abspath, join, normpath
import re

from stat import S_ISREG
from snakemake.remote import (
    AbstractRemoteObject,
    AbstractRemoteProvider,
    AbstractRemoteRetryObject,
)
from snakemake.exceptions import WorkflowError, XRootDFileException

try:
    from XRootD import client
    from XRootD.client.flags import DirListFlags, MkDirFlags, StatInfoFlags
except ImportError as e:
    raise WorkflowError(
        "The Python 3 package 'XRootD' must be installed to use XRootD "
        "remote() file functionality. %s" % e.msg
    )


class RemoteProvider(AbstractRemoteProvider):
    supports_default = True

    def __init__(
        self,
        *args,
        keep_local=False,
        stay_on_remote=False,
        is_default=False,
        url_decorator=lambda x: x,
        **kwargs,
    ):
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs,
        )

        self._xrd = XRootDHelper(url_decorator)
        self._url_decorator = url_decorator

    def remote_interface(self):
        return self._xrd

    def remote(self, *args, **kwargs):
        if "url_decorator" in kwargs:
            return super().remote(*args, **kwargs)
        else:
            return super().remote(*args, url_decorator=self._url_decorator, **kwargs)

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "root://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["root://", "roots://", "rootk://"]


class RemoteObject(AbstractRemoteRetryObject):
    """This is a class to interact with XRootD servers."""

    def __init__(
        self,
        *args,
        keep_local=False,
        stay_on_remote=False,
        provider=None,
        url_decorator=lambda x: x,
        **kwargs,
    ):
        super(RemoteObject, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            provider=provider,
            **kwargs,
        )

        if provider:
            self._xrd = provider.remote_interface()
        else:
            self._xrd = XRootDHelper()

        self._url_decorator = url_decorator

    def to_plainstr(self):
        return self._url_decorator(super().to_plainstr())

    # === Implementations of abstract class members ===

    def exists(self):
        return self._xrd.exists(self._url_decorator(self.remote_file()))

    def mtime(self):
        if self.exists():
            return self._xrd.file_last_modified(self._url_decorator(self.remote_file()))
        else:
            raise XRootDFileException(
                "The file does not seem to exist remotely: %s" % self.remote_file()
            )

    def size(self):
        if self.exists():
            return self._xrd.file_size(self._url_decorator(self.remote_file()))
        else:
            return self._iofile.size_local

    def _download(self):
        assert not self.stay_on_remote
        self._xrd.copy(self._url_decorator(self.remote_file()), self.file())

    def _upload(self):
        assert not self.stay_on_remote
        self._xrd.copy(self.file(), self._url_decorator(self.remote_file()))

    @property
    def name(self):
        return self.local_file()

    @property
    def list(self):
        dirname = os.path.dirname(self._iofile.constant_prefix()) + "/"
        files = list(self._xrd.list_directory_recursive(dirname))
        return [normpath(f) for f in files]

    def remove(self):
        self._xrd.remove(self.remote_file())


class XRootDHelper(object):
    def __init__(self, url_decorator=lambda x: x):
        self._clients = {}
        self._url_decorator = url_decorator

    def get_client(self, domain):
        try:
            return self._clients[domain]
        except KeyError:
            self._clients[domain] = client.FileSystem(domain)
            return self._clients[domain]

    def _parse_url(self, url):
        match = re.search(
            "(?P<domain>(?:[A-Za-z]+://)[A-Za-z0-9:@\_\-\.]+\:?/)(?P<path>.+)", url
        )
        if match is None:
            return None

        domain = match.group("domain")

        dirname, filename = os.path.split(match.group("path"))
        # We need a trailing / to keep XRootD happy
        dirname += "/"

        # and also make sure we supply a non-relative path
        # (snakemake removes double-slash // characters)
        if not dirname.startswith("/"):
            dirname = "/" + dirname

        return domain, dirname, filename

    def _parse_url_with_token(self, url):
        """Same as parse_url but also separates the token for the cases that need it"""
        ret = self._parse_url(url)
        if ret is None:
            # In the case we have a local file, not a ROOT URL
            dirname = os.path.dirname(url)
            filename = os.path.basename(url)
            return None, dirname, filename, None
        domain, dirname, filename = ret
        components = filename.rsplit("?", 1)
        filename = components[0]
        token = None
        if len(components) > 1:
            token = components[1]
        return domain, dirname, filename, token

    def exists(self, url):
        domain, dirname, filename = self._parse_url(url)

        status, statInfo = self.get_client(domain).stat(os.path.join(dirname, filename))

        if not status.ok:
            if status.errno == 3011:
                return False
            raise XRootDFileException(
                "Error stating URL "
                + os.path.join(dirname, filename)
                + " on domain "
                + domain
                + "\n"
                + repr(status)
                + "\n"
                + repr(statInfo)
            )

        return True
        # return not (
        #     (statInfo.flags & StatInfoFlags.IS_DIR)
        #     or (statInfo.flags & StatInfoFlags.OTHER)
        # )

    def _get_statinfo(self, url):
        domain, dirname, filename, token = self._parse_url_with_token(url)
        matches = [
            f for f in self.list_directory(domain, dirname) if f.name == filename
        ]

        assert len(matches) > 0
        if len(matches) > 1:
            # -- check matches for consistency
            # There is a transient effect in XRootD
            # where a file may match more than once.
            # This is okay as long as the statinfo
            # is the same for all of them.
            relevant_properties = [  # we only need to check front-facing attributes
                x
                for x in dir(matches[0].statinfo)
                if not (x[:1] == "_" or x[-2:] == "__")
            ]
            assert all(
                getattr(m.statinfo, p) == getattr(matches[0].statinfo, p)
                for m in matches[1:]
                for p in relevant_properties
            )

        return matches[0].statinfo

    def file_last_modified(self, filename):
        return self._get_statinfo(filename).modtime

    def file_size(self, filename):
        return self._get_statinfo(filename).size

    def copy(self, source, destination):
        # Prepare the source path for XRootD
        if not self._parse_url(source):
            source = abspath(source)
            filename = os.path.basename(source)
        else:
            domain, dirname, filename, token = self._parse_url_with_token(source)
            source = (
                f"{domain}/{dirname}/{filename}?{token}"
                if token
                else f"{domain}/{dirname}/{filename}"
            )

        # Prepare the destination path for XRootD
        if self._parse_url(destination):
            (
                dest_domain,
                dest_dirname,
                dest_filename,
                dest_token,
            ) = self._parse_url_with_token(destination)
            destination = (
                f"{dest_domain}/{dest_dirname}/{dest_filename}?{dest_token}"
                if dest_token
                else f"{dest_domain}/{dest_dirname}/{dest_filename}"
            )
            self.makedirs(dest_domain, dest_dirname)
        else:
            destination = abspath(destination)
            dest_filename = os.path.basename(destination)
            if not os.path.isdir(os.path.dirname(destination)):
                os.makedirs(os.path.dirname(destination))

        # Checking that we keep the same filename
        assert filename == dest_filename

        # Perform the copy operation
        process = client.CopyProcess()
        process.add_job(source, destination, force=True)
        process.prepare()
        status, returns = process.run()
        if not status.ok or not returns[0]["status"].ok:
            raise XRootDFileException(
                "Error copying from " + source + " to " + destination,
                repr(status),
                repr(returns),
            )

    def makedirs(self, domain, dirname):
        print("Making directories", domain, dirname)
        assert dirname.endswith("/")
        status, _ = self.get_client(domain).mkdir(
            self._url_decorator(dirname), MkDirFlags.MAKEPATH
        )
        if not status.ok:
            raise XRootDFileException(
                "Failed to create directory " + dirname, repr(status)
            )

    def list_directory(self, domain, dirname):
        status, dirlist = self.get_client(domain).dirlist(
            self._url_decorator(dirname), DirListFlags.STAT
        )
        if not status.ok:
            raise XRootDFileException(
                "Error listing directory "
                + dirname
                + " on domain "
                + domain
                + "\n"
                + repr(status)
                + "\n"
                + repr(dirlist)
            )
        return dirlist.dirlist

    def list_directory_recursive(self, start_url):
        assert start_url.endswith("/")
        domain, dirname, filename = self._parse_url(start_url)
        assert not filename
        filename = join(dirname, filename)
        for f in self.list_directory(domain, dirname):
            if f.statinfo.flags & StatInfoFlags.IS_DIR:
                for _f_name in self.list_directory_recursive(
                    domain + dirname + f.name + "/"
                ):
                    yield _f_name
            else:
                # Only yield files as directories don't have timestamps on XRootD
                yield domain + dirname + f.name

    def remove(self, url):
        domain, dirname, filename = self._parse_url(url)
        filename = join(dirname, filename)
        status, _ = self.get_client(domain).rm(self._url_decorator(filename))
        if not status.ok:
            raise XRootDFileException(
                "Failed to remove file "
                + filename
                + " from remote "
                + domain
                + "\n"
                + repr(status)
            )
