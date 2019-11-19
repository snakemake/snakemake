__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

import os
import requests
import re
import hashlib
from collections import namedtuple

from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError
from snakemake.common import lazy_property
from snakemake.utils import makedirs

ZenFileInfo = namedtuple("ZenFileInfo", ["checksum", "filesize", "id",  "download"])

class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, stay_on_remote=False, deposition=None, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, deposition=deposition, **kwargs)
        self._zen = ZENHelper(*args, deposition=deposition, **kwargs)
    
    def remote_interface(self):
        return self._zen

    @property
    def default_protocol(self):
        return "https://"
    
    @property
    def available_protocols(self):
        return ["http://", "https://"]

class RemoteObject(AbstractRemoteObject):
    def __init__(self, *args, keep_local=False, stay_on_remote=False, provider=None, deposition=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, stay_on_remote=stay_on_remote, provider=provider, deposition=deposition, **kwargs)
        if provider:
            self._zen = provider.remote_interface()
        else:
            self._zen = ZENHelper(*args, deposition=deposition, **kwargs)

    # === Implementations of abstract class members ===
    def _stats(self):
        return self._zen.get_files()[os.path.basename(self.local_file())]
    
    def exists(self):
        return os.path.basename(self.local_file()) in self._zen.get_files()

    def size(self):
        return self._stats().filesize

    def mtime(self):
        # There is no mtime info provided by Zenodo
        # Hence, the files are always considered to be "ancient".
        return 0

    def download(self):
        self._zen.download(self.local_file())

    def upload(self):
        self._zen.upload(self.local_file(), self.remote_file())

    @property
    def list(self):
        return [i for i in self._zen.get_files()]
    
    @property
    def name(self):
        return self.local_file()

class ZENHelper(object):
    def __init__(self, *args, deposition=None, **kwargs):
        if not deposition:
            raise IndexError("Deposition number must be supplied.")
        self.deposition = deposition

        if "access_token" in kwargs:
            self._access_token = kwargs.pop("access_token")
        else:
            try:
                self._access_token = os.environ["ZENODO_PAT"]
            except KeyError:
                raise WorkflowError("$ZENODO_PAT must be given as environment variable. "
                                    "Otherwise 'access_token' need to be passed in as kwargs.")

    def _api_request(self,
                    url,
                    method = "GET",
                    data = None,
                    headers = {},
                    files = None,
                    json = False):
        
        # Create a session with a hook to raise error on bad request.
        session = requests.Session()
        session.hooks = {"response": lambda r, *args, **kwargs: r.raise_for_status()}
        session.headers["Authorization"] = "Bearer {}".format(self._access_token)
        session.headers.update(headers)

        # Run query
        r = session.request(method =  method, url = url, data = data, files = files)
        if json:
            msg = r.json()
            return msg
        else:
            return r
    
    def get_files(self):
        headers = {"Content-Type": "application/json"}
        files = self._api_request(
            self._baseurl + "api/deposit/depositions/{}/files".format(self.deposition),
            headers = headers,
            json = True)
        return {os.path.basename(f["filename"]): 
            ZenFileInfo(f["checksum"], int(f["filesize"]), f["id"], f["links"]["download"]) 
            for f in files}

    def download(self, remote_file):
        # Get stats with download link
        stats = self.get_files()[os.path.basename(remote_file)]
        r = self._api_request(stats.download)
        
        local_md5 = hashlib.md5()

        # Make dir if missing
        makedirs(os.path.dirname(os.path.realpath(remote_file)))
        
        # Download file
        with open(remote_file, "wb") as rf:
            for chunk in r.iter_content(chunk_size = 1024 * 1024 * 10):
                local_md5.update(chunk)
                rf.write(chunk)
        local_md5 = local_md5.hexdigest()

        if local_md5 != stats.checksum:
            raise WorkflowError("File checksums do not match for remote file id: {}".format(stats.id))  
    
    def upload(self, local_file, remote_file):
        # Current stable API supports 100MB per file.
        with open(local_file, "rb") as lf:
            self._api_request(
                self._baseurl + "api/deposit/depositions/{}/files".format(self.deposition),
                method = "POST",
                data = {"filename": remote_file},
                files = {"file": lf})