"""Azure Blob Storage handling
"""

__author__ = "Sebastian Kurscheid"
__copyright__ = "Copyright 2021, Sebastian Kurscheid"
__email__ = "sebastian.kurscheid@anu.edu.au"
__license__ = "MIT"

# built-ins
import os
import re

# snakemake specific
# /

# module specific
from snakemake.exceptions import WorkflowError, AzureFileException
from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider

# service provider support
try:
    from azure.storage.blob import BlobServiceClient
    import azure.core.exceptions
except ImportError as e:
    raise WorkflowError(
        "The Python 3 package 'azure-storage-blob' "
        "need to be installed to use Azure Storage remote() file functionality. %s"
        % e.msg
    )


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

        self._as = AzureStorageHelper(*args, **kwargs)

    def remote_interface(self):
        return self._as

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "ab://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["ab://"]


class RemoteObject(AbstractRemoteObject):
    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(
            *args, keep_local=keep_local, provider=provider, **kwargs
        )

        if provider:
            self._as = provider.remote_interface()
        else:
            self._as = AzureStorageHelper(*args, **kwargs)

    # === Implementations of abstract class members ===
    def exists(self):
        if self._matched_as_path:
            return self._as.exists_in_container(self.container_name, self.blob_name)
        raise AzureFileException(
            "The file cannot be parsed as an Azure Blob path in form 'container/blob': %s"
            % self.local_file()
        )

    def mtime(self):
        if self.exists():
            # b = self.blob_service_client.get_blob_client(self.container_name, self.blob_name)
            # return b.get_blob_properties().last_modified
            t = self._as.blob_last_modified(self.container_name, self.blob_name)
            return t
        raise AzureFileException(
            "The file does not seem to exist remotely: %s" % self.local_file()
        )

    def size(self):
        if self.exists():
            return self._as.blob_size(self.container_name, self.blob_name)
        return self._iofile.size_local

    def download(self):
        if self.exists():
            os.makedirs(os.path.dirname(self.local_file()), exist_ok=True)
            self._as.download_from_azure_storage(
                self.container_name, self.blob_name, destination_path=self.local_file()
            )
            os.sync()
            return self.local_file()
        return None

    def upload(self):
        self._as.upload_to_azure_storage(
            container_name=self.container_name,
            blob_name=self.blob_name,
            file_path=self.local_file(),
        )

    @property
    def list(self):
        return self._as.list_blobs(self.container_name)

    # # === Related methods ===
    @property
    def _matched_as_path(self):
        return re.search(
            "(?P<container_name>[^/]*)/(?P<blob_name>.*)", self.local_file()
        )

    def as_create_stub(self):
        if self._matched_as_path:
            if not self.exists:
                self._as.download_from_azure_storage(
                    self.container_name,
                    self.blob_name,
                    self.file,
                    create_stub_only=True,
                )
        else:
            raise AzureFileException(
                "The file to be downloaded cannot be parsed as an Azure Storage path in form 'container/blob': %s"
                % self.local_file()
            )

    @property
    def container_name(self):
        "return container name component of the path"
        if len(self._matched_as_path.groups()) == 2:
            return self._matched_as_path.group("container_name")
        return None

    @property
    def name(self):
        return self.blob_name

    @property
    def blob_name(self):
        "return the blob name component of the path"
        if len(self._matched_as_path.groups()) == 2:
            return self._matched_as_path.group("blob_name")
        return None


# Actual Azure specific functions, adapted from S3.py


class AzureStorageHelper(object):
    def __init__(self, *args, **kwargs):
        if "stay_on_remote" in kwargs:
            del kwargs["stay_on_remote"]

        # if not handed down explicitely, try to read credentials from
        # environment variables.
        for (csavar, envvar) in [
            ("account_url", "AZ_BLOB_ACCOUNT_URL"),
            ("credential", "AZ_BLOB_CREDENTIAL"),
        ]:
            if csavar not in kwargs and envvar in os.environ:
                kwargs[csavar] = os.environ.get(envvar)
        assert (
            "account_url" in kwargs
        ), "Missing AZ_BLOB_ACCOUNT_URL env var (and possibly AZ_BLOB_CREDENTIAL)"
        # remove leading '?' from SAS if needed
        # if kwargs.get("sas_token", "").startswith("?"):
        #    kwargs["sas_token"] = kwargs["sas_token"][1:]

        # by right only account_key or sas_token should be set, but we let
        # BlobServiceClient deal with the ambiguity
        self.blob_service_client = BlobServiceClient(**kwargs)

    def container_exists(self, container_name):
        return any(
            True for _ in self.blob_service_client.list_containers(container_name)
        )

    def upload_to_azure_storage(
        self,
        container_name,
        file_path,
        blob_name=None,
        use_relative_path_for_blob_name=True,
        relative_start_dir=None,
        extra_args=None,
    ):
        """Upload a file to Azure Storage
        This function uploads a file to an Azure Storage Container as a blob.
        Args:
            container_name: the name of the Azure container to use
            file_path: The path to the file to upload.
            blob_name: The name to set for the blob on Azure. If not specified, this will default to the
                name of the file.
        Returns: The blob_name of the file on Azure if written, None otherwise
        """
        file_path = os.path.realpath(os.path.expanduser(file_path))

        assert container_name, "container_name must be specified"
        assert os.path.exists(file_path), (
            "The file path specified does not exist: %s" % file_path
        )
        assert os.path.isfile(file_path), (
            "The file path specified does not appear to be a file: %s" % file_path
        )

        container_client = self.blob_service_client.get_container_client(container_name)
        try:
            container_client.create_container()
        except azure.core.exceptions.ResourceExistsError:
            pass

        if not blob_name:
            if use_relative_path_for_blob_name:
                if relative_start_dir:
                    path_blob_name = os.path.relpath(file_path, relative_start_dir)
                else:
                    path_blob_name = os.path.relpath(file_path)
            else:
                path_blob_name = os.path.basename(file_path)
            blob_name = path_blob_name
        blob_client = container_client.get_blob_client(blob_name)

        # upload_blob fails, if blob exists
        if self.exists_in_container(container_name, blob_name):
            blob_client.delete_blob()
        try:
            with open(file_path, "rb") as data:
                blob_client.upload_blob(data, blob_type="BlockBlob")
            return blob_client.get_blob_properties().name
        except Exception as e:
            raise WorkflowError("Error in creating blob. %s" % str(e))
            # return None

    def download_from_azure_storage(
        self,
        container_name,
        blob_name,
        destination_path=None,
        expandBlobNameIntoDirs=True,
        make_dest_dirs=True,
        create_stub_only=False,
    ):
        """Download a file from Azure Storage
        This function downloads an object from a specified Azure Storage container.
        Args:
            container_name: the name of the Azure Storage container to use (container name only)
            destination_path: If specified, the file will be saved to this path, otherwise cwd.
            expandBlobNameIntoDirs: Since Azure blob names can include slashes, if this is True (defult)
                then Azure blob names with slashes are expanded into directories on the receiving end.
                If it is False, the blob name is passed to os.path.basename() to get the substring
                following the last slash.
            make_dest_dirs: If this is True (default) and the destination path includes directories
                that do not exist, they will be created.
        Returns:
            The destination path of the downloaded file on the receiving end, or None if the destination_path
            could not be downloaded
        """
        assert container_name, "container_name must be specified"
        assert blob_name, "blob_name must be specified"
        if destination_path:
            destination_path = os.path.realpath(os.path.expanduser(destination_path))
        else:
            if expandBlobNameIntoDirs:
                destination_path = os.path.join(os.getcwd(), blob_name)
            else:
                destination_path = os.path.join(
                    os.getcwd(), os.path.basename(blob_name)
                )
        # if the destination path does not exist
        if make_dest_dirs:
            os.makedirs(os.path.dirname(destination_path), exist_ok=True)
        b = self.blob_service_client.get_blob_client(container_name, blob_name)
        if not create_stub_only:
            with open(destination_path, "wb") as my_blob:
                blob_data = b.download_blob()
                blob_data.readinto(my_blob)
        else:
            # just create an empty file with the right timestamps
            ts = b.get_blob_properties().last_modified.timestamp()
            with open(destination_path, "wb") as fp:
                os.utime(fp.name, (ts, ts))
        return destination_path

    def delete_from_container(self, container_name, blob_name):
        """Delete a file from Azure Storage container

        This function deletes an object from a specified Azure Storage container.

        Args:
            container_name: the name of the Azure Storage container to use (container name only, not endpoint)
            blob_name: the name of the blob to delete from the container

        Returns:
            nothing
        """
        assert container_name, "container_name must be specified"
        assert blob_name, "blob_name must be specified"
        b = self.blob_service_client.get_blob_client(container_name, blob_name)
        b.delete_blob()

    def exists_in_container(self, container_name, blob_name):
        """Returns whether the blob exists in the container

        Args:
            container_name: the name of the Azure Storage container (container name only, not endpoint)
            blob_name: the blob_name of the object to delete from the container

        Returns:
            True | False
        """

        assert (
            container_name
        ), 'container_name must be specified (did you try to write to "root" or forgot to set --default-remote-prefix?)'
        assert blob_name, "blob_name must be specified"
        cc = self.blob_service_client.get_container_client(container_name)
        return any(True for _ in cc.list_blobs(name_starts_with=blob_name))

    def blob_size(self, container_name, blob_name):
        """Returns the size of a blob

        Args:
            container_name: the name of the Azure Storage container (container name only, not endpoint)
            blob_name: the blob_name of the object to delete from the container

        Returns:
            Size in kb
        """
        assert container_name, "container_name must be specified"
        assert blob_name, "blob_name must be specified"

        b = self.blob_service_client.get_blob_client(container_name, blob_name)
        return b.get_blob_properties().size // 1024

    def blob_last_modified(self, container_name, blob_name):
        """Returns a timestamp of a blob

        Args:
            container_name: the name of the Azure Storage container (container name only, not endpoint)
            blob_name: the blob_name of the object to delete from the container

        Returns:
            timestamp
        """
        assert container_name, "container_name must be specified"
        assert blob_name, "blob_name must be specified"
        b = self.blob_service_client.get_blob_client(container_name, blob_name)
        return b.get_blob_properties().last_modified.timestamp()

    def list_blobs(self, container_name):
        """Returns a list of blobs from the container

        Args:
            container_name: the name of the Azure Storage container (container name only, not endpoint)

        Returns:
            list of blobs
        """
        assert container_name, "container_name must be specified"
        c = self.blob_service_client.get_container_client(container_name)
        return [b.name for b in c.list_blobs()]
