from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union
import re


from snakemake_interface_common.exceptions import WorkflowError


PATHVAR_NAME_REGEX = re.compile(r"^[a-z][a-z0-9_]+$")
PATHVAR_REGEX = re.compile(r"<(?P<name>[a-z][a-z0-9_]+)>")


class Pathvars:

    def __init__(self) -> None:
        self.items: Dict[str, str]
        self.level: Dict[str, int]

    @classmethod
    def from_other(cls, other: "Pathvars") -> "Pathvars":
        instance = cls()
        instance.items = dict(other.items)
        instance.level = dict(other.level)
        return instance

    @classmethod
    def from_raw(cls, items: Dict[str, str], level: Optional[int] = None) -> "Pathvars":
        cls.check_dict(items)
        instance = cls()
        instance.items = items
        if level is None:
            assert not items
            instance.level = {}
        else:
            instance.level = {key: level for key in items}
        return instance

    @classmethod
    def with_defaults(cls) -> "Pathvars":
        items = {
            "results": "results",
            "resources": "resources",
            "logs": "logs",
            "benchmarks": "benchmarks",
        }
        return cls.from_raw(
            items=items,
            level=4,
        )

    @classmethod
    def from_config(cls, config: Union[Mapping[Any, Any], Sequence[Any]]) -> "Pathvars":
        if isinstance(config, dict):
            config_pathvars = config.get("pathvars")
            if config_pathvars:
                # ignore type here as checked above
                return cls.from_raw(items=config_pathvars, level=1)  # type: ignore[arg-type]
        return cls.from_raw(items={})

    @classmethod
    def from_module(cls, module_pathvars: Dict[str, str]) -> "Pathvars":
        return cls.from_raw(items=module_pathvars, level=2)

    @classmethod
    def from_rule(cls, rule_pathvars: Dict[str, str]) -> "Pathvars":
        return cls.from_raw(items=rule_pathvars, level=0)

    @classmethod
    def from_workflow(cls, workflow_pathvars: Dict[str, str]) -> "Pathvars":
        return cls.from_raw(items=workflow_pathvars, level=3)

    def update(self, other: "Pathvars") -> None:
        for key, value in other.items.items():
            if key not in self.items or other.level[key] <= self.level[key]:
                self.items[key] = value
                self.level[key] = other.level[key]

    def get(self, name: str) -> str:
        return self.items[name]

    def apply(self, path: str) -> str:
        while PATHVAR_REGEX.search(path):
            path = PATHVAR_REGEX.sub(lambda item: self.get(item.group("name")), path)
        return path

    @classmethod
    def check_dict(cls, items: Dict[Any, Any]) -> None:
        errors = {}
        for key, value in items.items():
            if (
                not isinstance(key, str)
                or not PATHVAR_NAME_REGEX.match(key)
                or not isinstance(value, str)
            ):
                errors[key] = value
        if errors:
            errors = ",".join(f"{key}:{value}" for key, value in errors.items())
            raise WorkflowError(
                f"Pathvars have to be a mapping of str to str, with keys being valid pathvar names (i.e. alphanumeric + _, lower-case, no leading number). The following entries are invalid: {errors}"
            )
