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
    A typed file wrapper that serializes/deserializes a typed object (Dataclass or NamedTuple) as JSON.
    """

    def __init__(self, file: str, type_: Type[T]):
        self.file = file
        self.type_ = type_

    def __str__(self):
        return str(self.file)

    @classmethod
    def factory(cls, type_: Type[T]):
        """
        Returns a constructor for TypedFile bound to the given type.
        type_ should be a Dataclass, NamedTuple, or any class with asdict().
        """
        return lambda file: cls(file, type_)

    def dump(self, *args, **kwargs):
        """Construct an instance of type_ and serialize it to the file as JSON."""
        obj = self.type_(*args, **kwargs)
        # obj1 = self.type_(**typed_to_dict(obj))  # is it needed to verify the type?
        with open(str(self.file), "w") as f:
            json.dump(typed_to_dict(obj), fp=f)

    def load(self) -> T:
        """Deserialize the file and return an instance of type_."""
        with open(str(self.file), "r") as f:
            data = json.load(f)
        return self.type_(**data)

    @property
    def value(self):
        """Read the current value from disk."""
        return self.load()

    def __getattr__(self, name: str):
        """
        Transparently proxy attribute access to the deserialized object.
        Allows input.meta.some_threshold without an explicit .load() call.
        Only triggered when the attribute is not found on TypedFile itself.
        """
        return getattr(self.load(), name)
