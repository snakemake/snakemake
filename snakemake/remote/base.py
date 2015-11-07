
import re

from snakemake.remote import AbstractRemoteObject

class DomainObject(AbstractRemoteObject):
    """This is a mixin related to parsing components
        out of a location path specified as
        (host|IP):port/remote/location
    """
    def __init__(self, *args, **kwargs):
            super(DomainObject, self).__init__(*args, **kwargs)

    @property
    def _matched_address(self):
        return re.search("^(?P<host>[A-Za-z0-9\-\.]+)(?:\:(?P<port>[0-9]+))?(?P<path_remainder>.*)$", self._iofile._file)

    @property
    def name(self):
        return self.path_remainder
    
    @property
    def protocol(self):
        if self._matched_address:
            return self._matched_address.group("protocol")

    @property
    def host(self):
        if self._matched_address:
            return self._matched_address.group("host")

    @property
    def port(self):
        return self._matched_address.group("port")
    
    @property
    def path_prefix(self):
        # this is the domain and port, however specified before the path remainder
        return self._iofile._file[:self._iofile._file.index(self.path_remainder)]
    
    @property
    def path_remainder(self):
        if self._matched_address:
            return self._matched_address.group("path_remainder")

    @property
    def local_path(self):
        return self._iofile._file

    @property
    def remote_path(self):
        return self.path_remainder