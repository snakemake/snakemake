from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import typing


class A(ABC):
    @property
    @abstractmethod
    def test(self) -> typing.Optional[int]:
        ...


@dataclass
class B(A):
    test: int = frozenset()


print(B(test=2))
print(isinstance(B(test=2), A))
