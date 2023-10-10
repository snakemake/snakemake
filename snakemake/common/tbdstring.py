# A string that prints as TBD
# whatever interaction happens on this class, <TBD> shall be returned
class TBDString(str):
    # the second arg is necessary to avoid problems when pickling
    def __new__(cls, _=None):
        return str.__new__(cls, "<TBD>")

    def __getitem__(self, __item):
        return self

    def __bool__(self):
        return False

    def __add__(self, __other):
        return self

    def __sub__(self, __other):
        return self

    def __mul__(self, __other):
        return self

    def __matmul__(self, __other):
        return self

    def __truediv__(self, __other):
        return self

    def __floordiv__(self, __other):
        return self

    def __mod__(self, __other):
        return self

    def __divmod__(self, __other):
        return self

    def __pow__(self, __other):
        return self

    def __lshift__(self, __other):
        return self

    def __rshift__(self, __other):
        return self

    def __and__(self, __other):
        return self

    def __xor__(self, __other):
        return self

    def __or__(self, __other):
        return self

    def __neg__(self):
        return self

    def __pos__(self):
        return self

    def __abs__(self):
        return self

    def __invert__(self):
        return self

    def __complex__(self):
        return self

    def __int__(self):
        return self

    def __float__(self):
        return self

    def __index__(self):
        return self

    def __round__(self, ndigits=0):
        return self

    def __trunc__(self):
        return self

    def __floor__(self):
        return self

    def __ceil__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, __exc_type, __exc_value, __traceback):
        return self
