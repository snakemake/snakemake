__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import os, sys
from contextlib import contextmanager
import pickle
import time
import threading

# third-party
import boto
from moto import mock_s3

# intra-module
from snakemake.remote.S3 import RemoteObject as S3RemoteObject, RemoteProvider as S3RemoteProvider
from snakemake.remote.S3 import S3Helper
from snakemake.decorators import dec_all_methods

def noop():
    pass

def pickled_moto_wrapper(func):
    """
        This is a class decorator that in turn decorates all methods within
        a class to mock out boto calls with moto-simulated ones.
        Since the moto backends are not presistent across calls by default, 
        the wrapper also pickles the bucket state after each function call,
        and restores it before execution. This way uploaded files are available
        for follow-on tasks. Since snakemake may execute with multiple threads
        it also waits for the pickled bucket state file to be available before
        loading it in. This is a hackey alternative to using proper locks,
        but works ok in practice.
    """
    def wrapper_func(self, *args, **kwargs):
        moto_context_file = "motoState.p"

        moto_context = mock_s3()

        # load moto buckets from pickle
        if os.path.isfile(moto_context_file) and os.path.getsize(moto_context_file) > 0:
            with file_lock(moto_context_file):
                with open( moto_context_file, "rb" ) as f:
                    moto_context.backends["global"].buckets = pickle.load( f )

        moto_context.backends["global"].reset = noop

        mocked_function = moto_context(func)

        retval = mocked_function(self, *args, **kwargs)

        with file_lock(moto_context_file):
            with open( moto_context_file, "wb" ) as f:
                pickle.dump(moto_context.backends["global"].buckets, f)

        return retval
    return wrapper_func

@dec_all_methods(pickled_moto_wrapper, prefix=None)
class RemoteProvider(S3RemoteProvider):
    def __init__(self, *args, **kwargs):
        return super(RemoteProvider, self).__init__(*args, **kwargs)
        
@dec_all_methods(pickled_moto_wrapper, prefix=None)
class RemoteObject(S3RemoteObject):
    """ 
        This is a derivative of the S3 remote provider that mocks
        out boto-based S3 calls using the "moto" Python package.
        Only the initializer is different; it "uploads" the input 
        test file to the moto-simulated bucket at the start.
    """

    def __init__(self, *args, keep_local=False, provider=None, **kwargs):
        bucket_name = 'test-remote-bucket'
        test_file = "test.txt"

        conn = boto.connect_s3()
        if bucket_name not in [b.name for b in conn.get_all_buckets()]:
            conn.create_bucket(bucket_name)

        # "Upload" files that should be in S3 before tests...
        s3c = S3Helper()
        if not s3c.exists_in_bucket(bucket_name, test_file):
            s3c.upload_to_s3(bucket_name, test_file)

        return super(RemoteObject, self).__init__(*args, keep_local=keep_local, provider=provider, **kwargs)


# ====== Helpers =====

@contextmanager
def file_lock(filepath):
    lock_file = filepath + ".lock"

    while os.path.isfile(lock_file):
        time.sleep(0.1)

    with open(lock_file, 'w') as f:
        f.write("1")

    try:
        yield
    finally:
        if os.path.isfile(lock_file):
            os.remove(lock_file)

