import textwrap
import shutil


class Linter:
    def __init__(self, workflow, items):
        self.items = items
        self.workflow = workflow

    def lint(self, json=False):
        for item in self.items:
            item_lints = [
                lint for lint_item in self.lints() for lint in lint_rule(rule)
            ]

            if json:
                json_lints = []
            if item_lints:
                if json:
                    json_lints.append(
                        {
                            "for": self.item_desc_json(item),
                            "lints": [lint.__dict__ for lint in item_lints],
                        }
                    )
                else:
                    logger.info(
                        "Lints for {}:\n{}\n".format(
                            self.item_desc_plain(item),
                            "\n".join(map("    * {}".format, item_lints)),
                        )
                    )

            if json:
                return json_lints

    @abstractmethod
    def item_desc_json(self, item):
        pass

    @abstractmethod
    def item_desc_plain(self, item):
        pass

    def lints(self):
        return (
            method for name, method in self.__dict__.items() if name.startswith("lint_")
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
