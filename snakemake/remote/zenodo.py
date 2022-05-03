__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

import os
import hashlib
from collections import namedtuple
import requests
from requests.exceptions import HTTPError
from snakemake.remote import (
    AbstractRemoteObject,
    AbstractRemoteProvider,
    AbstractRemoteRetryObject,
)
from snakemake.exceptions import ZenodoFileException, WorkflowError
from snakemake.common import lazy_property


ZenFileInfo = namedtuple("ZenFileInfo", ["checksum", "filesize", "id", "download"])


class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, stay_on_remote=False, **kwargs):
        super(RemoteProvider, self).__init__(
            *args, stay_on_remote=stay_on_remote, **kwargs
        )
        self._zen = ZENHelper(*args, **kwargs)

    def remote_interface(self):
        return self._zen

    @property
    def default_protocol(self):
        return "https://"

    @property
    def available_protocols(self):
        return ["https://"]


class RemoteObject(AbstractRemoteRetryObject):
    def __init__(
        self, *args, keep_local=False, stay_on_remote=False, provider=None, **kwargs
    ):
        super(RemoteObject, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            provider=provider,
            **kwargs,
        )
        if provider:
            self._zen = provider.remote_interface()
        else:
            self._zen = ZENHelper(*args, **kwargs)

    # === Implementations of abstract class members ===
    def _stats(self):
        return self._zen.get_files()[os.path.basename(self.local_file())]

    def exists(self):
        return os.path.basename(self.local_file()) in self._zen.get_files()

    def size(self):
        if self.exists():
            return self._stats().filesize
        else:
            return self._iofile.size_local

    def mtime(self):
        # There is no mtime info provided by Zenodo.
        # Hence, the files are always considered to be "ancient".
        return 0

    def _download(self):
        stats = self._stats()
        download_url = stats.download
        r = self._zen._api_request(download_url)

        local_md5 = hashlib.md5()

        # Download file.
        with open(self.local_file(), "wb") as rf:
            for chunk in r.iter_content(chunk_size=1024 * 1024 * 10):
                local_md5.update(chunk)
                rf.write(chunk)
        local_md5 = local_md5.hexdigest()

        if local_md5 != stats.checksum:
            raise ZenodoFileException(
                "File checksums do not match for remote file id: {}".format(stats.id)
            )

    def _upload(self):
        with open(self.local_file(), "rb") as lf:
            self._zen._api_request(
                self._zen.bucket + "/{}".format(os.path.basename(self.remote_file())),
                method="PUT",
                data=lf,
            )

    @property
    def list(self):
        return [i for i in self._zen.get_files()]

    @property
    def name(self):
        return self.local_file()


class ZENHelper(object):
    def __init__(self, *args, **kwargs):

        try:
            self._access_token = kwargs.pop("access_token")
        except KeyError:
            raise WorkflowError(
                "Zenodo personal access token must be passed in as 'access_token' argument.\n"
                "Separate registration and access token is needed for Zenodo sandbox "
                "environment at https://sandbox.zenodo.org."
            )

        self.restricted_access_token = None
        self._restricted_access_cookies = None

        if "restricted_access_token" in kwargs:
            self.restricted_access_token = kwargs["restricted_access_token"]

        if "sandbox" in kwargs:
            self._sandbox = kwargs.pop("sandbox")
        else:
            self._sandbox = False

        if self._sandbox:
            self._baseurl = "https://sandbox.zenodo.org"
        else:
            self._baseurl = "https://zenodo.org"

        if "deposition" in kwargs:
            self.deposition = kwargs.pop("deposition")
            self.bucket = self.get_bucket()
        else:
            # Creating a new deposition, as deposition id was not supplied.
            self.deposition, self.bucket = self.create_deposition().values()

    def _api_request(
        self, url, method="GET", data=None, headers={}, files=None, json=False
    ):

        # Create a session with a hook to raise error on bad request.
        session = requests.Session()
        session.hooks = {"response": lambda r, *args, **kwargs: r.raise_for_status()}
        session.headers["Authorization"] = "Bearer {}".format(self._access_token)
        session.headers.update(headers)

        cookies = self.restricted_access_cookies

        # Run query.
        try:
            r = session.request(
                method=method, url=url, data=data, files=files, cookies=cookies
            )
            if json:
                msg = r.json()
                return msg
            else:
                return r
        except HTTPError as e:
            raise WorkflowError("Failed to connect to zenodo", e)

    def create_deposition(self):
        resp = self._api_request(
            method="POST",
            url=self._baseurl + "/api/deposit/depositions",
            headers={"Content-Type": "application/json"},
            data="{}",
            json=True,
        )
        return {"id": resp["id"], "bucket": resp["links"]["bucket"]}

    def get_bucket(self):
        resp = self._api_request(
            self._baseurl + "/api/deposit/depositions/{}".format(self.deposition),
            headers={"Content-Type": "application/json"},
            json=True,
        )
        return resp["links"]["bucket"]

    def get_files(self):
        files = self._api_request(
            self._baseurl + "/api/deposit/depositions/{}/files".format(self.deposition),
            headers={"Content-Type": "application/json"},
            json=True,
        )
        return {
            os.path.basename(f["filename"]): ZenFileInfo(
                f["checksum"], int(f["filesize"]), f["id"], f["links"]["download"]
            )
            for f in files
        }

    @property
    def restricted_access_cookies(self):
        """Retrieve cookies necessary for restricted access.

        Inspired by https://gist.github.com/slint/d47fe5628916d14b8d0b987ac45aeb66
        """
        if self.restricted_access_token and self._restricted_access_cookies is None:
            url = (
                self._baseurl
                + f"/record/{self.deposition}?token={self.restricted_access_token}"
            )
            resp = self._api_request(url)
            if "session" in resp["cookies"]:
                self._restricted_access_cookies = resp["cookies"]
            else:
                raise WorkflowError(
                    "Failure to retrieve session cookie with given restricted access token. "
                    f"Is the token valid? Please check by opening {url} manually in your browser."
                )
        return self._restricted_access_cookies
