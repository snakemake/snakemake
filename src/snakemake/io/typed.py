import json
import dataclasses
from typing import Generic, Type, TypeVar

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


class TypedFile(Generic[T]):
    """
    A wrapper to indicate that a file has a specific type.
    TODO: is it necessary to zip the file?
    """

    def __init__(self, file: str, type_: Type[T]):
        self.file = file
        self.type_ = type_

    def __str__(self):
        return str(self.file)

    @classmethod
    def factory(cls, type_: Type[T]):
        """
        type_ should be Dataclass or NamedTuple, or any class that can be set as json and reloaded
        """
        return lambda file: cls(file, type_)

    def dump(self, *args, **kwargs):
        obj = self.type_(*args, **kwargs)
        # obj1 = self.type_(**typed_to_dict(obj))  # is it needed to verify the type?
        with open(str(self.file), "w") as f:
            json.dump(typed_to_dict(obj), fp=f)

    def load(self) -> T:
        """Guard that the value can only be read after it has been written."""
        with open(str(self.file), "r") as f:
            data = json.load(f)
        return self.type_(**data)

    @property
    def value(self):
        return self.load()
