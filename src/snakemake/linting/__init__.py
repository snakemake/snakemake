import textwrap
import shutil
import inspect
from abc import ABC, abstractmethod

from snakemake.logging import logger

# (?!\\+) is a negative lookahead, which removes trailing
# '+'s from the match. There is a minor risk, that a user
# intentionally uses file name (parts) with a trailing '+'
# intentionally. The regex extension _should_ allow simple
# regexes as '\s+' in place of a tab separator to pass.
NAME_PATTERN = "[a-zA-Z_][a-zA-Z_0-9]*(?!\\+)"


class Linter(ABC):
    def __init__(self, workflow, items):
        self.items = items
        self.workflow = workflow

    def read_item(self, item):
        return item

    def lint(self, json=False):
        json_lints = [] if json else None
        linted = False
        for item in self.items:
            item_lints = [
                lint
                for lint_item in self.lints()
                for lint in lint_item(self.read_item(item))
            ]
            if not item_lints:
                continue
            linted = True
            if json:
                json_lints.append(
                    {
                        "for": self.item_desc_json(item),
                        "lints": [lint.__dict__ for lint in item_lints],
                    }
                )
            else:
                logger.warning(
                    "Lints for {}:\n{}\n".format(
                        self.item_desc_plain(item),
                        "\n".join(map("    * {}".format, item_lints)),
                    )
                )
        return json_lints, linted

    @abstractmethod
    def item_desc_json(self, item):
        pass

    @abstractmethod
    def item_desc_plain(self, item):
        pass

    def lints(self):
        return (
            method
            for name, method in inspect.getmembers(self)
            if name.startswith("lint_")
        )


class Lint:
    def __init__(self, title, body, links=None):
        self.title = title
        self.body = body
        self.links = links or []

    def __str__(self) -> str:
        width, _ = shutil.get_terminal_size()
        output = "{}:\n{}".format(
            self.title,
            "\n".join(
                map("      {}".format, textwrap.wrap(self.body, max(width - 6, 20)))
            ),
        )
        if self.links:
            output += "\n      Also see:\n{}".format(
                "\n".join(map("      {}".format, self.links))
            )
        return output
