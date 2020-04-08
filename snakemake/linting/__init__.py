import textwrap
import shutil
import inspect
from abc import ABC, abstractmethod

from snakemake.logging import logger

NAME_PATTERN = "[a-zA-Z_][a-zA-Z_0-9]*"


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
    def __init__(self, title, body, links=[]):
        self.title = title
        self.body = body
        self.links = links

    def __str__(self):
        width, _ = shutil.get_terminal_size()
        return "{}:\n{}\n      Also see:\n{}".format(
            self.title,
            "\n".join(
                map("      {}".format, textwrap.wrap(self.body, max(width - 6, 20)))
            ),
            "\n".join(map("      {}".format, self.links)),
        )
