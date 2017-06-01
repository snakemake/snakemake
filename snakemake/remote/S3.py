__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import os
import re
import math
import email.utils
import functools
import concurrent.futures

# module-specific
from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError, S3FileException

try:
    # third-party modules
    import boto
    from boto.s3.key import Key
    from filechunkio import FileChunkIO
except ImportError as e:
    raise WorkflowError("The Python 3 packages 'boto' and 'filechunkio' " +
        "need to be installed to use S3 remote() file functionality. %s" % e.msg)


class RemoteProvider(AbstractRemoteProvider):
    def __init__(self, *args, stay_on_remote=False, **kwargs):
        super(RemoteProvider, self).__init__(*args, stay_on_remote=stay_on_remote, **kwargs)

        self._s3c = S3Helper(*args, **kwargs)

    def remote_interface(self):
        return self._s3c

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return 's3://'

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ['s3://']


class RemoteObject(AbstractRemoteObject):
    """ This is a class to interact with the AWS S3 object store.
    """

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)

        if provider:
            self._s3c = provider.remote_interface()
        else:
            self._s3c = S3Helper(*args, **kwargs)

    # === Implementations of abstract class members ===

    def exists(self):
        if self._matched_s3_path:
            return self._s3c.exists_in_bucket(self.s3_bucket, self.s3_key)
        else:
            raise S3FileException("The file cannot be parsed as an s3 path in form 'bucket/key': %s" % self.local_file())

    def mtime(self):
        if self.exists():
            return self._s3c.key_last_modified(self.s3_bucket, self.s3_key)
        else:
            raise S3FileException("The file does not seem to exist remotely: %s" % self.local_file())

    def size(self):
        if self.exists():
            return self._s3c.key_size(self.s3_bucket, self.s3_key)
        else:
            return self._iofile.size_local

    def download(self):
        self._s3c.download_from_s3(self.s3_bucket, self.s3_key, self.local_file())
        os.sync() # ensure flush to disk

    def upload(self):
        if self.size() > 10 * 1024 * 1024: # S3 complains if multipart uploads are <10MB
            self._s3c.upload_to_s3_multipart(self.s3_bucket, self.local_file(), self.s3_key, encrypt_key=self.kwargs.get("encrypt_key", None))
        else:
            self._s3c.upload_to_s3(self.s3_bucket, self.local_file(), self.s3_key, encrypt_key=self.kwargs.get("encrypt_key", None))

    @property
    def list(self):
        return [k.name for k in self._s3c.list_keys(self.s3_bucket)]

    # === Related methods ===

    @property
    def _matched_s3_path(self):
        return re.search("(?P<bucket>[^/]*)/(?P<key>.*)", self.local_file())

    @property
    def s3_bucket(self):
        if len(self._matched_s3_path.groups()) == 2:
            return self._matched_s3_path.group("bucket")
        return None

    @property
    def name(self):
        return self.s3_key

    @property
    def s3_key(self):
        if len(self._matched_s3_path.groups()) == 2:
            return self._matched_s3_path.group("key")

    def s3_create_stub(self):
        if self._matched_s3_path:
            if not self.exists:
                self._s3c.download_from_s3(self.s3_bucket, self.s3_key, self.file, create_stub_only=True)
        else:
            raise S3FileException("The file to be downloaded cannot be parsed as an s3 path in form 'bucket/key': %s" %
                                  self.local_file())


class S3Helper(object):

    def __init__(self, *args, **kwargs):
        # as per boto, expects the environment variables to be set:
        # AWS_ACCESS_KEY_ID
        # AWS_SECRET_ACCESS_KEY
        # Otherwise these values need to be passed in as kwargs

        # allow key_id and secret to be specified with aws_, gs_, or no prefix.
        # Standardize to the aws_ prefix expected by boto.
        if "gs_access_key_id" in kwargs:
            kwargs["aws_access_key_id"] = kwargs.pop("gs_access_key_id")
        if "gs_secret_access_key" in kwargs:
            kwargs["aws_secret_access_key"] = kwargs.pop("gs_secret_access_key")
        if "access_key_id" in kwargs:
            kwargs["aws_access_key_id"] = kwargs.pop("access_key_id")
        if "secret_access_key" in kwargs:
            kwargs["aws_secret_access_key"] = kwargs.pop("secret_access_key")

        self.conn = boto.connect_s3(*args, **kwargs)

    def upload_to_s3(
            self,
            bucket_name,
            file_path,
            key=None,
            use_relative_path_for_key=True,
            relative_start_dir=None,
            replace=False,
            reduced_redundancy=False,
            headers=None,
            encrypt_key=False):
        """ Upload a file to S3

            This function uploads a file to an AWS S3 bucket.

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                file_path: The path to the file to upload.
                key: The key to set for the file on S3. If not specified, this will default to the
                    name of the file.
                use_relative_path_for_key: If set to True (default), and key is None, the S3 key will include slashes
                    representing the path of the file relative to the CWD. If False only the
                    file basename will be used for the key.
                relative_start_dir: The start dir to use for use_relative_path_for_key. No effect if key is set.
                replace: If True a file with the same key will be replaced with the one being written
                reduced_redundancy: Sets the file to AWS reduced redundancy storage.
                headers: additional heads to pass to AWS

            Returns: The key of the file on S3 if written, None otherwise
        """
        file_path = os.path.realpath(os.path.expanduser(file_path))

        assert bucket_name, "bucket_name must be specified"
        assert os.path.exists(file_path), "The file path specified does not exist: %s" % file_path
        assert os.path.isfile(file_path), "The file path specified does not appear to be a file: %s" % file_path

        try:
            b = self.conn.get_bucket(bucket_name)
        except:
            b = self.conn.create_bucket(bucket_name)

        k = Key(b)

        if key:
            k.key = key
        else:
            if use_relative_path_for_key:
                if relative_start_dir:
                    path_key = os.path.relpath(file_path, relative_start_dir)
                else:
                    path_key = os.path.relpath(file_path)
            else:
                path_key = os.path.basename(file_path)
            k.key = path_key
        try:
            bytes_written = k.set_contents_from_filename(
                file_path,
                replace=replace,
                reduced_redundancy=reduced_redundancy,
                headers=headers,
                encrypt_key=encrypt_key)
            if bytes_written:
                return k.key
            else:
                return None
        except:
            return None

    def download_from_s3(
            self,
            bucket_name,
            key,
            destination_path=None,
            expandKeyIntoDirs=True,
            make_dest_dirs=True,
            headers=None, create_stub_only=False):
        """ Download a file from s3

            This function downloads an object from a specified AWS S3 bucket.

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                destination_path: If specified, the file will be saved to this path, otherwise cwd.
                expandKeyIntoDirs: Since S3 keys can include slashes, if this is True (defult)
                    then S3 keys with slashes are expanded into directories on the receiving end.
                    If it is False, the key is passed to os.path.basename() to get the substring
                    following the last slash.
                make_dest_dirs: If this is True (default) and the destination path includes directories
                    that do not exist, they will be created.
                headers: Additional headers to pass to AWS

            Returns:
                The destination path of the downloaded file on the receiving end, or None if the destination_path
                could not be downloaded
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucket_name)
        k = Key(b)

        if destination_path:
            destination_path = os.path.realpath(os.path.expanduser(destination_path))
        else:
            if expandKeyIntoDirs:
                destination_path = os.path.join(os.getcwd(), key)
            else:
                destination_path = os.path.join(os.getcwd(), os.path.basename(key))

        # if the destination path does not exist
        if make_dest_dirs:
            os.makedirs(os.path.dirname(destination_path), exist_ok=True)

        k.key = key if key else os.path.basename(destination_path)

        try:
            if not create_stub_only:
                k.get_contents_to_filename(destination_path, headers=headers)
            else:
                # just create an empty file with the right timestamps
                with open(destination_path, 'wb') as fp:
                    modified_tuple = email.utils.parsedate_tz(k.last_modified)
                    modified_stamp = int(email.utils.mktime_tz(modified_tuple))
                    os.utime(fp.name, (modified_stamp, modified_stamp))
            return destination_path
        except:
            return None

    def _upload_part(self, bucket_name, multipart_id, part_num, source_path, offset, bytes_to_write, number_of_retries=5):

        def _upload(retries_remaining=number_of_retries):
            try:
                b = self.conn.get_bucket(bucket_name)
                for mp in b.get_all_multipart_uploads():
                    if mp.id == multipart_id:
                        with FileChunkIO(source_path, 'r', offset=offset, bytes=bytes_to_write) as fp:
                            mp.upload_part_from_file(fp=fp, part_num=part_num)
                        break
            except Exception() as e:
                if retries_remaining:
                    _upload(retries_remaining=retries_remaining - 1)
                else:
                    raise e

        _upload()

    def upload_to_s3_multipart(
            self,
            bucket_name,
            file_path,
            key=None,
            use_relative_path_for_key=True,
            relative_start_dir=None,
            replace=False,
            reduced_redundancy=False,
            headers=None,
            parallel_processes=4,
            encrypt_key=False):
        """ Upload a file to S3

            This function uploads a file to an AWS S3 bucket.

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                file_path: The path to the file to upload.
                key: The key to set for the file on S3. If not specified, this will default to the
                    name of the file.
                use_relative_path_for_key: If set to True (default), and key is None, the S3 key will include slashes
                    representing the path of the file relative to the CWD. If False only the
                    file basename will be used for the key.
                relative_start_dir: The start dir to use for use_relative_path_for_key. No effect if key is set.
                replace: If True a file with the same key will be replaced with the one being written
                reduced_redundancy: Sets the file to AWS reduced redundancy storage.
                headers: additional heads to pass to AWS
                parallel_processes: Number of concurrent uploads

            Returns: The key of the file on S3 if written, None otherwise
        """
        file_path = os.path.realpath(os.path.expanduser(file_path))

        assert bucket_name, "bucket_name must be specified"
        assert os.path.exists(file_path), "The file path specified does not exist: %s" % file_path
        assert os.path.isfile(file_path), "The file path specified does not appear to be a file: %s" % file_path

        try:
            b = self.conn.get_bucket(bucket_name)
        except:
            b = self.conn.create_bucket(bucket_name)

        path_key = None
        if key:
            path_key = key
        else:
            if use_relative_path_for_key:
                if relative_start_dir:
                    path_key = os.path.relpath(file_path, relative_start_dir)
                else:
                    path_key = os.path.relpath(file_path)
            else:
                path_key = os.path.basename(file_path)

        mp = b.initiate_multipart_upload(path_key, headers=headers, encrypt_key=encrypt_key)

        source_size = os.stat(file_path).st_size

        bytes_per_chunk = 52428800  # 50MB = 50 * 1024 * 1024
        chunk_count = int(math.ceil(source_size / float(bytes_per_chunk)))

        with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_processes) as executor:
            for i in range(chunk_count):
                offset = i * bytes_per_chunk
                remaining_bytes = source_size - offset
                bytes_to_write = min([bytes_per_chunk, remaining_bytes])
                part_num = i + 1
                executor.submit(functools.partial(self._upload_part, bucket_name, mp.id, part_num, file_path, offset, bytes_to_write))

        if len(mp.get_all_parts()) == chunk_count:
            mp.complete_upload()
            try:
                key = b.get_key(path_key)
                return key.key
            except:
                return None
        else:
            mp.cancel_upload()
            return None

    def delete_from_bucket(self, bucket_name, key, headers=None):
        """ Delete a file from s3

            This function deletes an object from a specified AWS S3 bucket.

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket
                headers: Additional headers to pass to AWS

            Returns:
                The name of the object deleted
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucket_name)
        k = Key(b)
        k.key = key
        ret = k.delete(headers=headers)
        return ret.name

    def exists_in_bucket(self, bucket_name, key, headers=None):
        """ Returns whether the key exists in the bucket

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket
                headers: Additional headers to pass to AWS

            Returns:
                True | False
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucket_name)
        k = Key(b)
        k.key = key
        return k.exists(headers=headers)

    def key_size(self, bucket_name, key, headers=None):
        """ Returns the size of a key based on a HEAD request

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket
                headers: Additional headers to pass to AWS

            Returns:
                Size in kb
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucket_name)
        k = b.lookup(key)

        return k.size

    def key_last_modified(self, bucket_name, key, headers=None):
        """ Returns a timestamp of a key based on a HEAD request

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket
                headers: Additional headers to pass to AWS

            Returns:
                timestamp
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucket_name)
        k = b.lookup(key)

        # email.utils parsing of timestamp mirrors boto whereas
        # time.strptime() can have TZ issues due to DST
        modified_tuple = email.utils.parsedate_tz(k.last_modified)
        epochTime = int(email.utils.mktime_tz(modified_tuple))

        return epochTime

    def list_keys(self, bucket_name):
        return self.conn.get_bucket(bucket_name).list()
