from typing import FrozenSet, Set, TypeVar, Union


T = TypeVar("T")
AnySet = Union[Set[T], FrozenSet[T]]
