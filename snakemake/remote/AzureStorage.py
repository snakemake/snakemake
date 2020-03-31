__author__ = "Sebastian Kurscheid"
__copyright__ = "Copyright 2019, Sebastian Kurscheid"
__email__ = "sebastian.kurscheid@anu.edu.au"
__license__ = "MIT"

# built-ins
import os
import re
import math
import functools
import concurrent.futures

# snakemake specific
from snakemake.common import lazy_property

# module specific
from snakemake.exceptions import WorkflowError, AzureFileException
from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider

# service provider support
try:
    from azure.storage.common.cloudstorageaccount import (
        CloudStorageAccount as AzureStorageAccount,
    )
except ImportError as e:
    raise WorkflowError(
        "The Python 3 packages 'azure-storage' and 'azure-storage-common' "
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
        else:
            raise AzureFileException(
                "The file cannot be parsed as an Azure Blob path in form 'container/blob': %s"
                % self.local_file()
            )

    def mtime(self):
        if self.exists():
            t = self._as.blob_last_modified(self.container_name, self.blob_name)
            return t
        else:
            raise AzureFileException(
                "The file does not seem to exist remotely: %s" % self.local_file()
            )

    def size(self):
        if self.exists():
            return self._as.blob_size(self.container_name, self.blob_name)
        else:
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
        if len(self._matched_as_path.groups()) == 2:
            return self._matched_as_path.group("container_name")
        return None

    @property
    def name(self):
        return self.blob_name

    @property
    def blob_name(self):
        if len(self._matched_as_path.groups()) == 2:
            return self._matched_as_path.group("blob_name")


# Actual Azure specific functions, adapted from S3.py


class AzureStorageHelper(object):
    def __init__(self, *args, **kwargs):
        if "stay_on_remote" in kwargs:
            del kwargs["stay_on_remote"]

        self.azure = AzureStorageAccount(**kwargs).create_block_blob_service()

    def container_exists(self, container_name):
        try:
            self.azure.exists(container_name=container_name)
            return True
        except:
            return False

    def upload_to_azure_storage(
        self,
        container_name,
        file_path,
        blob_name=None,
        use_relative_path_for_blob_name=True,
        relative_start_dir=None,
        extra_args=None,
    ):
        """ Upload a file to Azure Storage
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

        if not self.azure.exists(container_name):
            self.azure.create_container(container_name=container_name)
        if not blob_name:
            if use_relative_path_for_blob_name:
                if relative_start_dir:
                    path_blob_name = os.path.relpath(file_path, relative_start_dir)
                else:
                    path_blob_name = os.path.relpath(file_path)
            else:
                path_blob_name = os.path.basename(file_path)
            blob_name = path_blob_name
        b = self.azure
        try:
            b.create_blob_from_path(
                container_name, file_path=file_path, blob_name=blob_name
            )
            return b.get_blob_properties(container_name, blob_name=blob_name).name
        except:
            raise WorkflowError("Error in creating blob. %s" % e.msg)
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
        """ Download a file from Azure Storage
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
        b = self.azure
        try:
            if not create_stub_only:
                b.get_blob_to_path(
                    container_name=container_name,
                    blob_name=blob_name,
                    file_path=destination_path,
                )
            else:
                # just create an empty file with the right timestamps
                with open(destination_path, "wb") as fp:
                    os.utime(
                        fp.name,
                        (b.last_modified.timestamp(), b.last_modified.timestamp()),
                    )
            return destination_path
        except:
            return None

    def delete_from_container(self, container_name, blob_name):
        """ Delete a file from Azure Storage container

            This function deletes an object from a specified Azure Storage container.

            Args:
                container_name: the name of the Azure Storage container to use (container name only, not endpoint)
                blob_name: the name of the blob to delete from the container

            Returns:
                nothing
        """
        assert container_name, "container_name must be specified"
        assert blob_name, "blob_name must be specified"
        b = self.azure
        b.delete_blob(container_name, blob_name)

    def exists_in_container(self, container_name, blob_name):
        """ Returns whether the blob exists in the container

            Args:
                container_name: the name of the Azure Storage container (container name only, not endpoint)
                blob_name: the blob_name of the object to delete from the container

            Returns:
                True | False
        """
        assert container_name, "container_name must be specified"
        assert blob_name, "blob_name must be specified"
        try:
            return self.azure.exists(container_name, blob_name)
        except:
            return None

    def blob_size(self, container_name, blob_name):
        """ Returns the size of a blob

            Args:
                container_name: the name of the Azure Storage container (container name only, not endpoint)
                blob_name: the blob_name of the object to delete from the container

            Returns:
                Size in kb
        """
        assert container_name, "container_name must be specified"
        assert blob_name, "blob_name must be specified"

        try:
            b = self.azure.get_blob_properties(container_name, blob_name)
            return b.properties.content_length // 1024
        except:
            print("blob or container do not exist")
            return None

    def blob_last_modified(self, container_name, blob_name):
        """ Returns a timestamp of a blob

            Args:
                container_name: the name of the Azure Storage container (container name only, not endpoint)
                blob_name: the blob_name of the object to delete from the container

            Returns:
                timestamp
        """
        assert container_name, "container_name must be specified"
        assert blob_name, "blob_name must be specified"
        try:
            b = self.azure.get_blob_properties(container_name, blob_name)
            return b.properties.last_modified.timestamp()
        except:
            print("blob or container do not exist")
            return None

    def list_blobs(self, container_name):
        """ Returns a list of blobs from the container

            Args:
                container_name: the name of the Azure Storage container (container name only, not endpoint)

            Returns:
                list of blobs
        """
        assert container_name, "container_name must be specified"
        try:
            b = self.azure.list_blobs(container_name)
            return [o.name for o in b]
        except:
            print("Did you provide a valid container_name?")
            return None
