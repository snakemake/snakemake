__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import os
import math
import time
import email.utils
from time import mktime
import datetime
from multiprocessing import Pool

# third-party modules
import boto
from boto.s3.key import Key
from filechunkio import FileChunkIO


class S3Helper(object):

    def __init__(self, *args, **kwargs):
        # as per boto, expects the environment variables to be set:
        # AWS_ACCESS_KEY_ID
        # AWS_SECRET_ACCESS_KEY
        # Otherwise these values need to be passed in as kwargs
        self.conn = boto.connect_s3(*args, **kwargs)

    def upload_to_s3(
            self,
            bucketName,
            filePath,
            key=None,
            useRelativePathForKey=True,
            relativeStartDir=None,
            replace=False,
            reduced_redundancy=False,
            headers=None):
        """ Upload a file to S3

            This function uploads a file to an AWS S3 bucket.

            Args:
                bucketName: the name of the S3 bucket to use (bucket name only, not ARN)
                filePath: The path to the file to upload.
                key: The key to set for the file on S3. If not specified, this will default to the
                    name of the file.
                useRelativePathForKey: If set to True (default), and key is None, the S3 key will include slashes
                    representing the path of the file relative to the CWD. If False only the
                    file basename will be used for the key.
                relativeStartDir: The start dir to use for useRelativePathForKey. No effect if key is set.
                replace: If True a file with the same key will be replaced with the one being written
                reduced_redundancy: Sets the file to AWS reduced redundancy storage.
                headers: additional heads to pass to AWS

            Returns: The key of the file on S3 if written, None otherwise
        """
        filePath = os.path.realpath(os.path.expanduser(filePath))

        assert bucketName, "bucketName must be specified"
        assert os.path.exists(filePath), "The file path specified does not exist: %s" % filePath
        assert os.path.isfile(filePath), "The file path specified does not appear to be a file: %s" % filePath

        try:
            b = self.conn.get_bucket(bucketName)
        except:
            b = self.conn.create_bucket(bucketName)

        k = Key(b)

        if key:
            k.key = key
        else:
            if useRelativePathForKey:
                if relativeStartDir:
                    pathKey = os.path.relpath(filePath, relativeStartDir)
                else:
                    pathKey = os.path.relpath(filePath)
            else:
                pathKey = os.path.basename(filePath)
            k.key = pathKey
        try:
            bytesWritten = k.set_contents_from_filename(
                filePath,
                replace=replace,
                reduced_redundancy=reduced_redundancy,
                headers=headers)
            if bytesWritten:
                return k.key
            else:
                return None
        except:
            return None

    def download_from_s3(
            self,
            bucketName,
            key,
            destinationPath=None,
            expandKeyIntoDirs=True,
            makeDestDirs=True,
            headers=None, createStubOnly=False):
        """ Download a file from s3

            This function downloads an object from a specified AWS S3 bucket.

            Args:
                bucketName: the name of the S3 bucket to use (bucket name only, not ARN)
                destinationPath: If specified, the file will be saved to this path, otherwise cwd.
                expandKeyIntoDirs: Since S3 keys can include slashes, if this is True (defult)
                    then S3 keys with slashes are expanded into directories on the receiving end.
                    If it is False, the key is passed to os.path.basename() to get the substring
                    following the last slash.
                makeDestDirs: If this is True (default) and the destination path includes directories
                    that do not exist, they will be created.
                headers: Additional headers to pass to AWS

            Returns:
                The destination path of the downloaded file on the receiving end, or None if the filePath
                could not be downloaded
        """
        assert bucketName, "bucketName must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucketName)
        k = Key(b)

        if destinationPath:
            destinationPath = os.path.realpath(os.path.expanduser(destinationPath))
        else:
            if expandKeyIntoDirs:
                destinationPath = os.path.join(os.getcwd(), key)
            else:
                destinationPath = os.path.join(os.getcwd(), os.path.basename(key))

        # if the destination path does not exist
        if not os.path.exists(os.path.dirname(destinationPath)) and makeDestDirs:
            os.makedirs(os.path.dirname(destinationPath))

        k.key = key if key else os.path.basename(filePath)

        try:
            if not createStubOnly:
                k.get_contents_to_filename(destinationPath, headers=headers)
            else:
                # just create an empty file with the right timestamps
                with open(destinationPath, 'wb') as fp:
                    modified_tuple = email.utils.parsedate_tz(k.last_modified)
                    modified_stamp = int(email.utils.mktime_tz(modified_tuple))
                    os.utime(fp.name, (modified_stamp, modified_stamp))
            return destinationPath
        except:
            return None

    def _upload_part(self, bucketName, multipart_id, part_num, source_path, offset, bytesToWrite, numberOfRetries=5):

        def _upload(retriesRemaining=numberOfRetries):
            try:
                b = self.conn.get_bucket(bucketName)
                for mp in b.get_all_multipart_uploads():
                    if mp.id == multipart_id:
                        with FileChunkIO(source_path, 'r', offset=offset, bytes=bytesToWrite) as fp:
                            mp.upload_part_from_file(fp=fp, part_num=part_num)
                        break
            except Exception() as e:
                if retriesRemaining:
                    _upload(retriesRemaining=retriesRemaining - 1)
                else:
                    raise e

        _upload()

    def upload_to_s3_multipart(
            self,
            bucketName,
            filePath,
            key=None,
            useRelativePathForKey=True,
            relativeStartDir=None,
            replace=False,
            reduced_redundancy=False,
            headers=None,
            parallel_processes=4):
        """ Upload a file to S3

            This function uploads a file to an AWS S3 bucket.

            Args:
                bucketName: the name of the S3 bucket to use (bucket name only, not ARN)
                filePath: The path to the file to upload.
                key: The key to set for the file on S3. If not specified, this will default to the
                    name of the file.
                useRelativePathForKey: If set to True (default), and key is None, the S3 key will include slashes
                    representing the path of the file relative to the CWD. If False only the
                    file basename will be used for the key.
                relativeStartDir: The start dir to use for useRelativePathForKey. No effect if key is set.
                replace: If True a file with the same key will be replaced with the one being written
                reduced_redundancy: Sets the file to AWS reduced redundancy storage.
                headers: additional heads to pass to AWS
                parallel_processes: Number of concurrent uploads

            Returns: The key of the file on S3 if written, None otherwise
        """
        filePath = os.path.realpath(os.path.expanduser(filePath))

        assert bucketName, "bucketName must be specified"
        assert os.path.exists(filePath), "The file path specified does not exist: %s" % filePath
        assert os.path.isfile(filePath), "The file path specified does not appear to be a file: %s" % filePath

        try:
            b = self.conn.get_bucket(bucketName)
        except:
            b = self.conn.create_bucket(bucketName)

        pathKey = None
        if key:
            pathKey = key
        else:
            if useRelativePathForKey:
                if relativeStartDir:
                    pathKey = os.path.relpath(filePath, relativeStartDir)
                else:
                    pathKey = os.path.relpath(filePath)
            else:
                pathKey = os.path.basename(filePath)

        mp = b.initiate_multipart_upload(pathKey, headers=headers)

        sourceSize = os.stat(filePath).st_size

        bytesPerChunk = 52428800  # 50MB = 50 * 1024 * 1024
        chunkCount = int(math.ceil(sourceSize / float(bytesPerChunk)))

        pool = Pool(processes=parallel_processes)
        for i in range(chunkCount):
            offset = i * bytesPerChunk
            remainingBytes = sourceSize - offset
            bytesToWrite = min([bytesPerChunk, remainingBytes])
            partNum = i + 1
            pool.apply_async(self._upload_part, [bucketName, mp.id, partNum, filePath, offset, bytesToWrite])
        pool.close()
        pool.join()

        if len(mp.get_all_parts()) == chunkCount:
            mp.complete_upload()
            try:
                key = b.get_key(pathKey)
                return key.key
            except:
                return None
        else:
            mp.cancel_upload()
            return None

    def delete_from_bucket(self, bucketName, key, headers=None):
        """ Delete a file from s3

            This function deletes an object from a specified AWS S3 bucket.

            Args:
                bucketName: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket
                headers: Additional headers to pass to AWS

            Returns:
                The name of the object deleted
        """
        assert bucketName, "bucketName must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucketName)
        k = Key(b)
        k.key = key
        ret = k.delete(headers=headers)
        return ret.name

    def exists_in_bucket(self, bucketName, key, headers=None):
        """ Returns whether the key exists in the bucket

            Args:
                bucketName: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket
                headers: Additional headers to pass to AWS

            Returns:
                True | False
        """
        assert bucketName, "bucketName must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucketName)
        k = Key(b)
        k.key = key
        return k.exists(headers=headers)

    def key_size(self, bucketName, key, headers=None):
        """ Returns the size of a key based on a HEAD request

            Args:
                bucketName: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket
                headers: Additional headers to pass to AWS

            Returns:
                Size in kb
        """
        assert bucketName, "bucketName must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucketName)
        k = b.lookup(key)

        return k.size

    def key_last_modified(self, bucketName, key, headers=None):
        """ Returns a timestamp of a key based on a HEAD request

            Args:
                bucketName: the name of the S3 bucket to use (bucket name only, not ARN)
                key: the key of the object to delete from the bucket
                headers: Additional headers to pass to AWS

            Returns:
                timestamp
        """
        assert bucketName, "bucketName must be specified"
        assert key, "Key must be specified"

        b = self.conn.get_bucket(bucketName)
        k = b.lookup(key)

        # email.utils parsing of timestamp mirrors boto whereas
        # time.strptime() can have TZ issues due to DST
        modified_tuple = email.utils.parsedate_tz(k.last_modified)
        epochTime = int(email.utils.mktime_tz(modified_tuple))

        return epochTime

    def list_keys(self, bucketName):
        return self.conn.get_bucket(bucketName).list()
