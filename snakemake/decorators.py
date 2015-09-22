__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import functools
import inspect


def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]

    return memoizer


def decAllMethods(decorator, prefix='test_'):

    def decClass(cls):
        for name, m in inspect.getmembers(cls, inspect.isfunction):
            if prefix == None or name.startswith(prefix):
                setattr(cls, name, decorator(m))
        return cls

    return decClass
