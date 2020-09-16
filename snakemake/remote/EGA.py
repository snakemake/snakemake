__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018-2019, Johannes Köster"
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


EGAFileInfo = namedtuple("EGAFileInfo", ["size", "status", "id", "checksum"])
EGAFile = namedtuple("EGAFile", ["dataset", "path"])


class RemoteProvider(AbstractRemoteProvider):
    def __init__(
        self,
        *args,
        keep_local=False,
        stay_on_remote=False,
        is_default=False,
        retry=5,
        **kwargs
    ):
        super().__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs
        )
        self.retry = retry
        self._token = None
        self._expires = None
        self._file_cache = dict()

    def _login(self):
        if self._expires is not None and self._expires > time.time():
            # token is still valid
            return

        # token will expire in 10 minutes
        # (we stop using it 10 seconds earlier to be sure)
        self._expires = time.time() + 10 * 60 * 60 - 10

        data = {
            "grant_type": "password",
            "client_id": self._client_id(),
            "scope": "openid",
            "client_secret": self._client_secret(),
            "username": self._username(),
            "password": self._password(),
        }
        for i in range(3):
            try:
                r = requests.post(
                    "https://ega.ebi.ac.uk:8443/ega-openid-connect-server/token",
                    headers={"Content-Type": "application/x-www-form-urlencoded"},
                    data=data,
                )
            except requests.exceptions.ConnectionError as e:
                time.sleep(5)
                if i == 2:
                    raise WorkflowError("Error contacting EGA.", e)

        if r.status_code != 200:
            raise WorkflowError("Login to EGA failed with:\n{}".format(r.text))
        r = r.json()
        # store session token
        try:
            self._token = r["access_token"]
        except KeyError:
            raise WorkflowError("Login to EGA failed:\n{}".format(r))

    def _expire_token(self):
        self._expires = None

    @property
    def token(self):
        self._login()
        return self._token

    def api_request(
        self,
        url_suffix,
        url_prefix="https://ega.ebi.ac.uk:8051/elixir/",
        json=True,
        post=False,
        **params
    ):
        """Make an API request.

        Args:
            url_suffix (str): Part of REST API URL right of https://ega.ebi.ac.uk:8051/elixir/
            params (dict): Parameters to pass, except session
        """

        url = url_prefix + url_suffix
        headers = (
            {"Accept": "application/json"}
            if json
            else {"Accept": "application/octet-stream"}
        )
        headers["Authorization"] = "Bearer {}".format(self.token)

        for i in range(3):
            try:
                if post:
                    r = requests.post(
                        url,
                        stream=not json,
                        data=params,
                        params={"session": self.token},
                        headers=headers,
                    )
                else:
                    params = dict(params)
                    params["session"] = self.token
                    r = requests.get(
                        url, stream=not json, params=params, headers=headers
                    )
            except requests.exceptions.ConnectionError as e:
                time.sleep(5)
                if i == 2:
                    raise WorkflowError("Error contacting EGA.", e)
        if r.status_code != 200:
            raise WorkflowError(
                "Access to EGA API endpoint {} failed " "with:\n{}".format(url, r.text)
            )
        if json:
            msg = r.json()
            return msg
        else:
            return r

    def get_files(self, dataset):
        if dataset not in self._file_cache:
            files = self.api_request(
                "data/metadata/datasets/{dataset}/files".format(dataset=dataset)
            )
            self._file_cache[dataset] = {
                os.path.basename(f["fileName"])[:-4]: EGAFileInfo(
                    int(f["fileSize"]), f["fileStatus"], f["fileId"], f["checksum"]
                )
                for f in files
            }
        return self._file_cache[dataset]

    @property
    def default_protocol(self):
        return "ega://"

    @property
    def available_protocols(self):
        return ["ega://"]

    @classmethod
    def _username(cls):
        return cls._credentials("EGA_USERNAME")

    @classmethod
    def _password(cls):
        return cls._credentials("EGA_PASSWORD")

    @classmethod
    def _client_id(cls):
        return cls._credentials("EGA_CLIENT_ID")

    @classmethod
    def _client_secret(cls):
        return cls._credentials("EGA_CLIENT_SECRET")

    @classmethod
    def _credentials(cls, name):
        try:
            return os.environ[name]
        except KeyError:
            raise WorkflowError(
                "$EGA_USERNAME, $EGA_PASSWORD, $EGA_CLIENT_ID, "
                "$EGA_CLIENT_SECRET must be given "
                "as environment variables."
            )


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
        stats = self._stats()

        r = self.provider.api_request(
            "data/files/{}?destinationFormat=plain".format(stats.id), json=False
        )

        local_md5 = hashlib.md5()

        # download file in chunks and calculate md5 on the fly
        os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)

        with open(self.local_file(), "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024 * 10):
                local_md5.update(chunk)
                f.write(chunk)
        local_md5 = local_md5.hexdigest()

        if local_md5 != stats.checksum:
            raise WorkflowError(
                "File checksums do not match for: {}".format(self.remote_file())
            )

    @lazy_property
    def parts(self):
        parts = self.local_file().split("/")
        if parts[0] != "ega":
            raise WorkflowError(
                "Invalid EGA remote file name. Must be 'ega/<dataset>/<filepath>'"
            )
        _, dataset, path = self.local_file().split("/", 2)
        return EGAFile(dataset, path)
