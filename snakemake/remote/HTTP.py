__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import os
import re
import collections
import email.utils
from contextlib import contextmanager

# module-specific
from snakemake.remote import AbstractRemoteProvider, DomainObject
from snakemake.exceptions import HTTPFileException, WorkflowError

try:
    # third-party modules
    import requests
except ImportError as e:
    raise WorkflowError("The Python 3 package 'requests' " +
        "must be installed to use HTTP(S) remote() file functionality. %s" % e.msg)


class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, stay_on_remote=False, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return 'https://'

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ['http://', 'https://']

    def remote(self, value, *args, insecure=None, **kwargs):
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, collections.Iterable):
            values = value
        else:
            raise TypeError('Invalid type ({}) passed to remote: {}'.format(type(value), value))

        for i, file in enumerate(values):
            match = re.match('^(https?)://.+', file)
            if match:
                protocol, = match.groups()
                if protocol == 'https' and insecure:
                    raise SyntaxError('insecure=True cannot be used with a https:// url')
                if protocol == 'http' and insecure not in [None, False]:
                    raise SyntaxError('insecure=False cannot be used with a http:// url')
            else:
                if insecure:
                    values[i] = 'http://' + file
                else:
                    values[i] = 'https://' + file

        return super(RemoteProvider, self).remote(values, *args, **kwargs)


class RemoteObject(DomainObject):
    """ This is a class to interact with an HTTP server.
    """

    def __init__(self, *args, keep_local=False, provider=None, additional_request_string="", **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)
        self.additional_request_string = additional_request_string

    # === Implementations of abstract class members ===

    @contextmanager #makes this a context manager. after 'yield' is __exit__()
    def httpr(self, verb="GET", stream=False):
        # if args have been provided to remote(), use them over those given to RemoteProvider()
        args_to_use = self.provider.args
        if len(self.args):
            args_to_use = self.args

        # use kwargs passed in to remote() to override those given to the RemoteProvider()
        # default to the host and port given as part of the file, falling back to one specified
        # as a kwarg to remote() or the RemoteProvider (overriding the latter with the former if both)
        kwargs_to_use = {}
        kwargs_to_use["username"] = None
        kwargs_to_use["password"] = None
        kwargs_to_use["auth"] = None

        for k,v in self.provider.kwargs.items():
            kwargs_to_use[k] = v
        for k,v in self.kwargs.items():
            kwargs_to_use[k] = v

        # Check that in case authentication kwargs are provided, they are either ("username", "password") combination
        # or "auth", but not both.
        if kwargs_to_use["username"] and kwargs_to_use["password"] and kwargs_to_use["auth"]:
            raise TypeError("Authentication accepts either username and password or requests.auth object")
        # If "username" and "password" kwargs are provided, use those to construct a tuple for "auth". Neither
        # requests.head() nor requests.get() accept them as-is.
        if kwargs_to_use["username"] and kwargs_to_use["password"]:
            kwargs_to_use["auth"] = (kwargs_to_use["username"], kwargs_to_use["password"])
        # Delete "username" and "password" from kwargs
        del kwargs_to_use["username"]
        del kwargs_to_use["password"]

        url = self.remote_file() + self.additional_request_string

        if verb.upper() == "GET":
            r = requests.get(url, *args_to_use, stream=stream, **kwargs_to_use)
        if verb.upper() == "HEAD":
            r = requests.head(url, *args_to_use, **kwargs_to_use)

        yield r
        r.close()

    def exists(self):
        if self._matched_address:
            with self.httpr(verb="HEAD") as httpr:
                # if a file redirect was found
                if httpr.status_code in range(300,308):
                    raise HTTPFileException("The file specified appears to have been moved (HTTP %s), check the URL or try adding 'allow_redirects=True' to the remote() file object: %s" % (httpr.status_code, httpr.url))
                return httpr.status_code == requests.codes.ok
            return False
        else:
            raise HTTPFileException("The file cannot be parsed as an HTTP path in form 'host:port/abs/path/to/file': %s" % self.local_file())

    def mtime(self):
        if self.exists():
            with self.httpr(verb="HEAD") as httpr:

                file_mtime = self.get_header_item(httpr, "last-modified", default=0)

                modified_tuple = email.utils.parsedate_tz(file_mtime)
                epochTime = email.utils.mktime_tz(modified_tuple)

                return epochTime
        else:
            raise HTTPFileException("The file does not seem to exist remotely: %s" % self.remote_file())

    def size(self):
        if self.exists():
            with self.httpr(verb="HEAD") as httpr:

                content_size = int(self.get_header_item(httpr, "content-size", default=0))

                return content_size
        else:
            return self._iofile.size_local

    def download(self, make_dest_dirs=True):
        with self.httpr(stream=True) as httpr:
            if self.exists():
                # if the destination path does not exist
                if make_dest_dirs:
                    os.makedirs(os.path.dirname(self.local_path), exist_ok=True)

                    with open(self.local_path, 'wb') as f:
                        for chunk in httpr.iter_content(chunk_size=1024):
                            if chunk: # filter out keep-alives
                                f.write(chunk)
                    os.sync() # ensure flush to disk
            else:
                raise HTTPFileException("The file does not seem to exist remotely: %s" % self.remote_file())

    def upload(self):
        raise HTTPFileException("Upload is not permitted for the HTTP remote provider. Is an output set to HTTP.remote()?")

    def get_header_item(self, httpr, header_name, default):
        """
            Since HTTP header capitalization may differ, this returns
            a header value regardless of case
        """

        header_value = default
        for k,v in httpr.headers.items():
            if k.lower() == header_name:
                header_value = v
        return header_value

    @property
    def list(self):
        raise HTTPFileException("The HTTP Remote Provider does not currently support list-based operations like glob_wildcards().")
