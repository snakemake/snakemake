from itertools import chain
import shutil
import textwrap

from snakemake.logging import logger
from snakemake.io import get_wildcard_names


def lint_rules(rules):
    for rule in rules:
        rule_lints = [lint for lint_rule in lints for lint in lint_rule(rule)]

        if rule_lints:
            logger.info(
                "Lints for rule {} (line {}, {}):\n{}\n".format(
                    rule.name,
                    rule.lineno,
                    rule.snakefile,
                    "\n".join(map("    * {}".format, rule_lints)),
                )
            )


def lint_params_prefix(rule):
    for param, value in rule.params.items():
        if isinstance(value, str) and any(
            f.startswith(value) for f in chain(rule.input, rule.output)
        ):
            yield Lint(
                title="Param {} is a prefix of input or output file but hardcoded".format(
                    param
                ),
                body="If this is meant to represent a file path prefix, it will fail when running "
                "workflow in environments without a shared filesystem. "
                "Instead, provide a function that infers the appropriate prefix from the input or "
                "output file, e.g.: lambda w, input: os.path.splitext(input[0])[0]",
                links=[
                    "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#non-file-parameters-for-rules",
                    "https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#tutorial-input-functions",
                ],
            )


def lint_log_directive(rule):
    if not rule.log and not rule.norun:
        yield Lint(
            title="No log directive defined",
            body="Without a log directive, all output will be printed "
            "to the terminal. In distributed environments, this means "
            "that errors are harder to discover. In local environments, "
            "output of concurrent jobs will be mixed and become unreadable.",
            links=[
                "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#log-files"
            ],
        )


def lint_not_used_params(rule, valid_names={"input", "output", "log", "params"}):
    if rule.shellcmd:
        for name in get_wildcard_names(rule.shellcmd):
            if name not in valid_names:
                yield Lint(
                    title="Shell command directly uses discouraged variable.",
                    body="It is recommended to pass all files as input and output, and non-file parameters "
                    "via the params directive. Otherwise, provenance tracking is less accurate.",
                    links=[
                        "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#non-file-parameters-for-rules"
                    ],
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


lints = [
    func
    for name, func in globals().items()
    if name.startswith("lint_") and not name == "lint_rules"
]

print(lints)
