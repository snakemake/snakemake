__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2021, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

import functools
import inspect


def dec_all_methods(decorator, prefix="test_"):
    def dec_class(cls):
        for name, m in inspect.getmembers(cls, inspect.isfunction):
            if prefix == None or name.startswith(prefix):
                setattr(cls, name, decorator(m))
        return cls

    return dec_class
