from typing import TypeVar, Union

T = TypeVar("T")
AnySet = Union[set[T], frozenset[T]]
