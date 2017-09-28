__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import os
import re
import math
import functools
from concurrent import futures
import contextlib
from io import BytesIO
import itertools
import tempfile
import gzip
import zlib
import shutil
from functools import partial

# module-specific
from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError, S3FileException
from snakemake.utils import available_cpu_count

try:
    # third-party modules
    import boto3
    import botocore
except ImportError as e:
    raise WorkflowError("The Python 3 package 'boto3' "
        "needs to be installed to use S3 remote() file functionality. %s" % e.msg)


def range_finder(length, range_size):
    if range_size>length:
        return (0,length)

    return [(lambda x: (x,(x+range_size)-1) if x+range_size<=length else (x,length))(x) for x in range(0,length+1)[::range_size]]

def grouper(src, size):
    src_iter = iter(src)
    while True:
        cur_chunk = list(itertools.islice(src_iter, size))
        if not cur_chunk:
            break
        yield cur_chunk

def grouped_ranges(length, chunk_size, chunks_per_group):
    return grouper(range_finder(length, chunk_size), chunks_per_group)


def ranged_get(bucket_name, key, kwargs, range):
    """
    Params:
        k: boto3 object
        range: two-tuple of the byte range to return (inclusive)
    """
    boto3_session = boto3.session.Session()
    s3 = boto3_session.resource('s3', **kwargs)
    obj = s3.Object(bucket_name, key)

    range_str = "bytes={}-{}".format(*range)
    res = obj.get(Range=range_str)["Body"].read()
    return res

class RemoteProvider(AbstractRemoteProvider):

    supports_default = True

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
        self._s3c.download_from_s3(self.s3_bucket, self.s3_key, self.local_file(), auto_decompress=self.kwargs.get("auto_decompress", True), streaming_decompress=self.kwargs.get("streaming_decompress", False))
        os.sync() # ensure flush to disk

    def upload(self):
        self._s3c.upload_to_s3(self.s3_bucket, self.local_file(), self.s3_key, extra_args=self.kwargs.get("ExtraArgs", None), config=self.kwargs.get("Config", None))

    @property
    def list(self):
        return self._s3c.list_keys(self.s3_bucket)

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

        self.init_kwargs = kwargs

        self.s3 = boto3.resource('s3', **kwargs)

    def bucket_exists(self, bucket_name):
        try:
            self.s3.meta.client.head_bucket(Bucket=bucket_name)
            return True
        except:
            return False

    def upload_to_s3(
            self,
            bucket_name,
            file_path,
            key=None,
            use_relative_path_for_key=True,
            relative_start_dir=None,
            extra_args=None,
            gzip_compress=False,
            config=None):
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

            Returns: The key of the file on S3 if written, None otherwise
        """
        file_path = os.path.realpath(os.path.expanduser(file_path))

        assert bucket_name, "bucket_name must be specified"
        assert os.path.exists(file_path), "The file path specified does not exist: %s" % file_path
        assert os.path.isfile(file_path), "The file path specified does not appear to be a file: %s" % file_path

        if not self.bucket_exists(bucket_name):
            self.s3.create_bucket(Bucket=bucket_name)

        if not key:
            if use_relative_path_for_key:
                if relative_start_dir:
                    path_key = os.path.relpath(file_path, relative_start_dir)
                else:
                    path_key = os.path.relpath(file_path)
            else:
                path_key = os.path.basename(file_path)
            key = path_key

        k = self.s3.Object(bucket_name, key)

        try:
            if gzip_compress:
                if not extra_args:
                    extra_args = {'ContentEncoding': 'gzip'}
                else:
                    extra_args['ContentEncoding'] = 'gzip'

                if not compressed_fp:
                    compressed_fp = BytesIO()
                with gzip.GzipFile(fileobj=compressed_fp, mode='wb') as gz:
                    shutil.copyfileobj(file_path, gz)
                compressed_fp.seek(0)
                k.upload_fileobj(compressed_fp, ExtraArgs=extra_args, Config=config)
            else:
                k.upload_file(file_path, ExtraArgs=extra_args, Config=config)
        except:
            raise
            return None

    def download_from_s3(
            self,
            bucket_name,
            key,
            destination_path=None,
            expandKeyIntoDirs=True,
            make_dest_dirs=True,
            create_stub_only=False,
            auto_decompress=True,
            streaming_decompress=False):
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

            Returns:
                The destination path of the downloaded file on the receiving end, or None if the destination_path
                could not be downloaded
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        b = self.s3.Bucket(bucket_name)

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

        k = self.s3.Object(bucket_name, key)
        try:
            if not create_stub_only:
                if auto_decompress:
                    if k.content_encoding == "gzip":  
                        if not streaming_decompress:                       
                            def streamingbody_iter(streamingbody, chunk_size=1024):
                                '''
                                    Construct a generator to return reasonably large data chunks
                                    from boto3's StreamingBody object, which by default
                                    returns everything (loading it into mem).
                                '''
                                chunk = streamingbody.read(amt=chunk_size)
                                while chunk[0] != b"": #bytearray():
                                    yield chunk
                                    chunk = streamingbody.read(amt=chunk_size)

                            obj_ranged_get = partial(ranged_get, bucket_name, key, self.init_kwargs)

                            with tempfile.TemporaryFile() as temp_out:
                                count = 0
                                decompressor = zlib.decompressobj(32 + zlib.MAX_WBITS)  # offset 32 to skip the header
                                with futures.ProcessPoolExecutor(max_workers=available_cpu_count()*2) as executor:
                                    for range_group in grouped_ranges(k.content_length, 20*1024*1024, int(available_cpu_count()*1.25)):
                                        for res in executor.map(obj_ranged_get, range_group, chunksize=1):
                                            decompressed_chunk = decompressor.decompress(res)
                                            count+=1
                                            if decompressed_chunk:
                                                temp_out.write(decompressed_chunk)
                                        os.sync()
                                        temp_out.flush()
                                shutil.move(temp_out, destination_path)                        
                        else:
                            # streaming decompression without runaway memory usage
                            # slower than chunked downloading
                            def stream_gzip_decompress(stream):
                                decompressor = zlib.decompressobj(32 + zlib.MAX_WBITS)  # offset 32 to skip the header
                                for chunk in stream:
                                    decompressed_chunk = decompressor.decompress(chunk)
                                    if decompressed_chunk:
                                        yield decompressed_chunk

                            def streamingbody_iter(streamingbody, chunk_size=20*1024*1024):
                                '''
                                    Construct a generator to return reasonably large data chunks
                                    from boto3's StreamingBody object, which by default
                                    returns everything (loading it into mem).
                                '''
                                chunk = streamingbody.read(amt=chunk_size)
                                while chunk[0] != bytearray():
                                    yield chunk
                                    chunk = streamingbody.read(amt=chunk_size)

                            fd,temp_out = tempfile.mkstemp()
                            count = 0
                            with open(temp_out,"wb") as outf:
                                for data in stream_gzip_decompress( streamingbody_iter(k.get()["Body"]) ):
                                    outf.write(data)
                                    outf.flush()
                            shutil.move(temp_out, destination_path)
                else:
                    k.download_file(destination_path)
            else:
                # just create an empty file with the right timestamps
                with open(destination_path, 'wb') as fp:
                    os.utime(fp.name, (k.last_modified.timestamp(), k.last_modified.timestamp()))
            return destination_path
        except:
            return None

    def delete_from_bucket(self, bucket_name, key):
        """ Delete a file from s3

            This function deletes an object from a specified AWS S3 bucket.

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket

            Returns:
                The name of the object deleted
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        k = self.s3.Object(bucket_name, key)
        ret = k.delete()
        return ret.name

    def exists_in_bucket(self, bucket_name, key):
        """ Returns whether the key exists in the bucket

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket

            Returns:
                True | False
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        try:
            self.s3.Object(bucket_name, key).load()
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] == "404":
                return False
            else:
                raise
        return True

    def key_size(self, bucket_name, key):
        """ Returns the size of a key based on a HEAD request

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket

            Returns:
                Size in kb
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        k = self.s3.Object(bucket_name, key)

        return k.content_length // 1024

    def key_last_modified(self, bucket_name, key):
        """ Returns a timestamp of a key based on a HEAD request

            Args:
                bucket_name: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket

            Returns:
                timestamp
        """
        assert bucket_name, "bucket_name must be specified"
        assert key, "Key must be specified"

        k = self.s3.Object(bucket_name, key)

        return k.last_modified.timestamp()

    def list_keys(self, bucket_name):
        b = self.s3.Bucket(bucket_name)
        return [o.key for o in b.objects]
