from dataclasses import dataclass, field
from typing import Any, Dict, Optional
import re


from snakemake_interface_common.exceptions import WorkflowError


PATHVAR_NAME_REGEX = re.compile(r"^[a-z][a-z0-9_]+$")
PATHVAR_REGEX = re.compile(r"<(?P<name>[a-z][a-z0-9_]+)>")


class Pathvars:

    def __init__(self, parent: Optional["Pathvars"] = None) -> None:
        self.items: Dict[str, str] = {}
        if not parent:
            # set defaults if nothing inherited
            self.items["results"] = "results"
            self.items["resources"] = "resources"
            self.items["logs"] = "logs"
            self.items["benchmarks"] = "benchmarks"
        else:
            self.items.update(parent.items)

    def register_config(self, config: Dict[Any, Any]) -> None:
        config_pathvars = config.get("pathvars")
        if config_pathvars:
            if not isinstance(config_pathvars, dict) or not all(isinstance(key, str) and PATHVAR_NAME_REGEX.match(key) and isinstance(value, str) for key, value in config_pathvars.items()):
                raise WorkflowError("The pathvars key in config files has to contain a mapping of str to str, with keys being valid pathvar names (i.e. alphanumeric + _, lower-case, no leading number)")
            self.items.update(config_pathvars)

    def get(self, name: str) -> str:
        try:
            return self.items[name]
        except KeyError:
            raise WorkflowError(f"Undefined pathvar {name}.")

    def register(self, **items: str) -> None:
        self.items.update(items)

    def spawn_new_namespace(self) -> "Pathvars":
        return Pathvars(parent=self)

    def apply(self, path: str) -> str:
        while PATHVAR_REGEX.search(path):
            path = PATHVAR_REGEX.sub(lambda item: self.get(item.group("name")), path)
        return path