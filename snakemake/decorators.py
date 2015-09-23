__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2015, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import functools
import inspect

def decAllMethods(decorator, prefix='test_'):

    def decClass(cls):
        for name, m in inspect.getmembers(cls, inspect.isfunction):
            if prefix == None or name.startswith(prefix):
                setattr(cls, name, decorator(m))
        return cls

    return decClass
