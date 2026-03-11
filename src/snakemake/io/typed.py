import dataclasses
import json
from typing import Callable, Generic, Type, TypeVar

T = TypeVar("T")


def typed_to_dict(obj):
    """
    Convert NamedTuple/Dataclass/[objects with asdict function] to json.
    All the key must be valid varnames and all the values must be json serializable.
    """
    if dataclasses.is_dataclass(obj):
        if not isinstance(obj, type):
            return dataclasses.asdict(obj)
    elif isinstance(obj, tuple):
        if hasattr(obj, "_asdict"):
            return obj._asdict()  # type: ignore[reportAttributeAccessIssue]
    elif hasattr(obj, "asdict"):
        return obj.asdict()  # type: ignore[reportAttributeAccessIssue]
    raise NotImplementedError(f"Cannot convert {type(obj)} to dict")


class TypedFile(str, Generic[T]):
    """
    A file path associated with a specific type for structured JSON serialization.

    Behaves as a plain string (file path) in all contexts, but provides ``dump()`` and ``load()`` methods to serialize and deserialize instances of the associated type to and from disk as JSON.

    The associated type must be convertible to a dict (see :func:`typed_to_dict`).
    Type annotations serve as documentation only; no runtime type checking is performed.
    """

    type_: Type[T]

    def __new__(cls, value: str, type_: Type[T]):
        self = super().__new__(cls, value)
        self.type_ = type_
        return self

    def dump(self, *args, **kwargs):
        """Construct an instance of type_ and serialize it to the file as JSON."""
        obj = self.type_(*args, **kwargs)
        # obj1 = self.type_(**typed_to_dict(obj))  # is it needed to verify the type?
        with open(str(self), "w") as f:
            json.dump(typed_to_dict(obj), fp=f)

    def load(self) -> T:
        """Deserialize the file and return an instance of type_."""
        with open(str(self), "r") as f:
            data = json.load(f)
        return self.type_(**data)


def typed_factory(type_: Type[T]) -> Callable[[str], TypedFile[T]]:
    """
    Returns a constructor for TypedFile bound to the given type.
    type_ should be a Dataclass, NamedTuple, or any class with asdict().
    """
    return lambda file: TypedFile(file, type_)
