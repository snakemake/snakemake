from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Union
import re


from snakemake_interface_common.exceptions import WorkflowError

from snakemake.exceptions import UndefinedPathvarException


PATHVAR_NAME_REGEX = re.compile(r"^[a-z][a-z0-9_]*$")
PATHVAR_REGEX = re.compile(r"<(?P<name>[a-z][a-z0-9_]*)>")


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
    def _from_raw(
        cls, items: Dict[Any, Any], level: Optional[int] = None
    ) -> "Pathvars":
        cls.check_dict(items)
        instance = cls()
        instance.items = items
        if level is None:
            assert not items
            instance.level = {}
        else:
            instance.level = {key: level for key in items}
        instance.check()
        return instance

    def check(self) -> None:
        # Ensure that there are no cyclic references as they would lead to infinite loops
        # during pathvar expansion.
        def dfs(key: str, seen: Set[str]) -> None:
            if key in seen:
                cycle = ", ".join(seen | {key})
                raise WorkflowError(
                    "Cyclic pathvar reference detected between: "
                    f"{cycle} (pathvars: {fmt_pathvars(self.items)})"
                )
            seen.add(key)
            value = self.items[key]
            for match in PATHVAR_REGEX.finditer(value):
                ref_key = match.group("name")
                if ref_key in self.items:
                    dfs(ref_key, seen)
            seen.remove(key)

        for key in self.items:
            dfs(key, set())

    @classmethod
    def with_defaults(cls) -> "Pathvars":
        items = {
            "results": "results",
            "resources": "resources",
            "logs": "logs",
            "benchmarks": "benchmarks",
        }
        return cls._from_raw(
            items=items,
            level=4,
        )

    @classmethod
    def from_config(
        cls, config: Union[Mapping[Any, Any], Sequence[Any]], module_level: bool = False
    ) -> "Pathvars":
        if isinstance(config, dict):
            config_pathvars = config.get("pathvars")
            if config_pathvars:
                # ignore type here as checked above
                return cls._from_raw(items=config_pathvars, level=1 if module_level else 2)  # type: ignore[arg-type]
        return cls._from_raw(items={})

    @classmethod
    def from_module(cls, module_pathvars: Dict[str, str]) -> "Pathvars":
        return cls._from_raw(items=module_pathvars, level=1)

    @classmethod
    def from_rule(cls, rule_pathvars: Dict[str, str]) -> "Pathvars":
        return cls._from_raw(items=rule_pathvars, level=0)

    @classmethod
    def from_workflow(cls, workflow_pathvars: Dict[str, str]) -> "Pathvars":
        return cls._from_raw(items=workflow_pathvars, level=3)

    def update(self, other: "Pathvars") -> None:
        for key, value in other.items.items():
            if key not in self.items or other.level[key] <= self.level[key]:
                self.items[key] = value
                self.level[key] = other.level[key]

        # ensure that the resulting pathvars are still valid after updating
        self.check()

    def get(self, name: str) -> str:
        return self.items[name]

    def apply(self, path: str) -> str:
        applied_path = path
        while PATHVAR_REGEX.search(applied_path):
            try:
                applied_path = PATHVAR_REGEX.sub(
                    lambda item: self.get(item.group("name")), applied_path
                )
            except KeyError as e:
                intermediate_msg = (
                    f" (intermediate path: {applied_path})."
                    if applied_path != path
                    else "."
                )
                raise UndefinedPathvarException(
                    f"Undefined pathvar '{e.args[0]}' when expanding "
                    f"pathvars in {path}{intermediate_msg}"
                ) from e
        return applied_path

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
            raise WorkflowError(
                "Pathvars have to be a mapping of str to str, with keys being "
                "valid pathvar names (i.e. alphanumeric + _, lower-case, no "
                "leading number). The following entries are "
                f"invalid: {fmt_pathvars(errors)}"
            )


def fmt_pathvars(items: Dict[str, str]) -> str:
    return ",".join(f"{key}='{value}'" for key, value in items.items())
