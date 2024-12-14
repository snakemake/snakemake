from typing import Union
from enum import Enum


def rule_item_factory(name: str):
    """Allows to access input, output etc. from statements inside
    a rule but outside of run/shell etc. blocks. Returns an object that
    defers evaluation to the DAG phase.
    """
    if name == "threads":

        def inner(_wildcards, threads):
            return threads

        return inner
    return RuleItemProxy(name)


class KeyKind(Enum):
    ATTRIBUTE = 0
    ITEM = 1


class RuleItemProxy:
    """Proxy class for deferring access to attributes and keys of the given item to DAG
    phase.
    """

    def __init__(self, name):
        self.name = name

    def __getattr__(self, name: str):
        return self._deferred_get(name, KeyKind.ATTRIBUTE)

    def __getitem__(self, key):
        return self._deferred_get(key, KeyKind.ITEM)

    def _deferred_get(self, key: Union[str, int], kind: KeyKind):
        if kind == KeyKind.ATTRIBUTE:

            def _get(item):
                return getattr(item, key)

        elif kind == KeyKind.ITEM:

            def _get(item):
                return item[key]

        else:
            raise ValueError("kind must be 'attribute' or 'item'")

        if self.name == "output":

            def inner(wildcards, output):
                return _get(output)

        elif self.name == "input":

            def inner(wildcards, input):
                return _get(input)

        elif self.name == "params":

            def inner(wildcards, params):
                return _get(params)

        elif self.name == "resources":

            def inner(wildcards, resources):
                return _get(resources)

        else:
            raise ValueError(f"Unknown item type: {self.name}")

        return inner
