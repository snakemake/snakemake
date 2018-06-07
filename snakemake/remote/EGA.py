__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018, Johannes Köster"
__email__ = "johannes.koester@tu-dortmund.de"
__license__ = "MIT"

import os
import json
import time
import uuid
from collections import namedtuple
import hashlib

import requests
from requests.auth import HTTPBasicAuth


from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError
from snakemake.common import lazy_property


try:
    import cryptography
except ImportError:
    raise WorkflowError("EGA remote provider needs cryptography to be installed.")


EGAFileInfo = namedtuple("EGAFileInfo", ["size", "status", "id"])
EGAFile = namedtuple("EGAFile", ["dataset", "path"])


class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, stay_on_remote=False, retry=5, **kwargs):
        super().__init__(*args, stay_on_remote=stay_on_remote, **kwargs)
        self.retry = retry
        self._token = None
        self._expires = None
        self._key = str(uuid.uuid4())[:32]
        self._file_cache = dict()

    def _login(self):
        if self._expires is not None and self._expires > time.time():
            # token is still valid
            return

        # token will expire in 10 minutes
        # (we stop using it 10 seconds earlier to be sure)
        self._expires = time.time() + 10 * 60 * 60 - 10
        auth = HTTPBasicAuth(self._username(), self._password())
        r = requests.get(
            "https://ega.ebi.ac.uk/ega/rest/access/v2/users/login",
            headers={"Accept": "application/json"},
            auth=auth)
        if r.status_code != 200:
            raise WorkflowError(
                "Login to EGA failed with:\n{}".format(r.text))
        r = r.json()
        # store session token
        result = r["response"]["result"]
        if result[0] == "success":
            self._token = result[1]
        else:
            raise WorkflowError(
                "Login to EGA failed with: {}".format(r["header"]["userMessage"]))

    def _expire_token(self):
        self._expires = None

    @property
    def token(self):
        self._login()
        return self._token

    def api_request(self,
                    url_suffix,
                    url_prefix="https://ega.ebi.ac.uk/ega/rest/access/v2/",
                    json=True,
                    post=False,
                    **params):
        """Make an API request.

        Args:
            url_suffix (str): Part of REST API URL right of https://ega.ebi.ac.uk/ega/rest/access/v2/
            params (dict): Parameters to pass, except session
        """

        url = url_prefix + url_suffix
        headers = {"Accept": "application/json"} if json else {
                   "Accept": "application/octet-stream"}

        for i in range(3):
            try:
                if post:
                    r = requests.post(url, stream=not json, data=params,
                                      params={"session": self.token},
                                      headers=headers)
                else:
                    params = dict(params)
                    params["session"] = self.token
                    r = requests.get(url, stream=not json,
                                     params=params, headers=headers)
            except requests.exceptions.ConnectionError as e:
                time.sleep(5)
                if i == 2:
                    raise WorkflowError("Error contacting EGA.", e)
        if r.status_code != 200:
            raise WorkflowError("Access to EGA API endpoint {} failed "
                                "with:\n{}".format(url, r.text))
        if json:
            msg = r.json()
            try:
                if msg["header"]["code"] == "991":
                    # session lost re-login
                    self._expire_token()
                    return self.api_request(
                        url_suffix,
                        url_prefix,
                        json=json,
                        post=post,
                        **params)
                return msg["response"]["result"]
            except KeyError as e:
                raise WorkflowError("Invalid response from EGA:\n{}".format(r.text))
        else:
            return r


    def get_files(self, dataset):
        if dataset not in self._file_cache:
            files = self.api_request(
                "datasets/{dataset}/files".format(dataset=dataset))
            self._file_cache[dataset] = {
                os.path.basename(f["fileName"])[:-4]:
                EGAFileInfo(int(f["fileSize"]), f["fileStatus"], f["fileID"])
                for f in files}
        return self._file_cache[dataset]

    @property
    def default_protocol(self):
        return "ega://"

    @property
    def available_protocols(self):
        return ["ega://"]

    @classmethod
    def _username(cls):
        try:
            return os.environ["EGA_USERNAME"]
        except KeyError:
            raise WorkflowError("$EGA_USERNAME and $EGA_PASSWORD must be given "
                                "as environment variables.")

    @classmethod
    def _password(cls):
        try:
            return os.environ["EGA_PASSWORD"]
        except KeyError:
            raise WorkflowError("$EGA_USERNAME and $EGA_PASSWORD must be given "
                                "as environment variables.")


class RemoteObject(AbstractRemoteObject):

    # === Implementations of abstract class members ===
    def _stats(self):
        return self.provider.get_files(self.parts.dataset)[self.parts.path]

    def exists(self):
        return self.parts.path in self.provider.get_files(self.parts.dataset)

    def size(self):
        return self._stats().size

    def mtime(self):
        # There is no mtime info provided by EGA
        # Hence, the files are always considered to be "ancient".
        return 0

    def download(self):
        # request file
        reid = str(uuid.uuid4())
        # create request
        r = self.provider.api_request(
            "requests/new/files/{}".format(self._stats().id),
            post=True,
            downloadrequest=json.dumps({"rekey": self.provider._key,
                             "downloadType": "STREAM",
                             "descriptor": reid}))

        try:
            # get ticket
            r = self.provider.api_request("requests/{}".format(reid))
            ticket = r[0]["ticket"]

            # download ticket (try 3 times)
            url_prefix = "http://ega.ebi.ac.uk/ega/rest/ds/v2/"
            for i in range(3):
                try:
                    r = self.provider.api_request("downloads/{}".format(ticket),
                                                  url_prefix=url_prefix,
                                                  json=False)
                except WorkflowError as e:
                    if i < 2:
                        continue
                    else:
                        raise e

            local_md5 = hashlib.md5()

            from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
            from cryptography.hazmat.backends import default_backend

            key = self.provider._key.encode()
            cipher = Cipher(algorithms.AES(key),
                            modes.CTR(key), backend=default_backend())
            decryptor = self.provider.cipher.decryptor()


            # download file in chunks, decrypt and calculate md5 on the fly
            os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)
            with open(self.local_file(), "wb") as f:
                for chunk in r.iter_content(chunk_size=512):
                    print(chunk)
                    local_md5.update(chunk)
                    f.write(decryptor.update(chunk))
                f.write(decryptor.finalize())
            local_md5 = local_md5.hexdigest()

            # check md5
            remote_md5, size = self.provider.api_request(
                "results/{}".format(ticket),
                url_prefix=url_prefix)

            if local_md5 != remote_md5:
                raise WorkflowError("File checksums do not match for: {}".format(
                    self.remote_file()))
        finally:
            self.provider.api_request("requests/delete/{}".format(reid))


    @lazy_property
    def parts(self):
        parts = self.local_file().split("/")
        if parts[0] != "ega":
            raise WorkflowError("Invalid EGA remote file name. Must be 'ega/<dataset>/<filepath>'")
        _, dataset, path = self.local_file().split("/", 2)
        return EGAFile(dataset, path)
