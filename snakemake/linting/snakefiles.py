import re
from itertools import chain

from snakemake.linting import Linter, Lint, links, NAME_PATTERN


class SnakefileLinter(Linter):
    def item_desc_plain(self, snakefile):
        return "snakefile {}".format(snakefile)

    def item_desc_json(self, snakefile):
        return {"snakefile": snakefile}
    
    def read_item(self, snakefile):
        return open(snakefile).read()

    def lint_absolute_paths(
        self, snakefile, regex=re.compile("(?P<quote>['\"])(?P<path>(?:/[^/]+?)+?)(?P=quote)")
    ):
        for match in regex.finditer(snakefile):
            line = get_line(match, snakefile)
            print(line)
            yield Lint(
                title="Absolute path \"{}\" in line {}".format(match.group("path"), line),
                body="Do not define absolute paths inside of the workflow, since this "
                "renders your workflow irreproducible on other machines. "
                "Use path relative to the working directory instead, or make the path "
                "configurable via a config file.",
                links=[links.config],
            )

    def lint_mixed_func_and_rules(
        self,
        snakefile,
        rule_regex=re.compile("rule .+?:"),
        func_regex=re.compile("def .+?:"),
    ):
        if rule_regex.search(snakefile) and func_regex.search(snakefile):
            yield Lint(
                title="Mixed rules and functions in same snakefile.",
                body="Small one-liner functions used only once should be "
                "defined as lambda expressions. Other functions should be collected "
                "in a common module, e.g. 'rules/common.smk'. This makes the workflow "
                "steps more readable.",
                links=[links.includes],
            )

    def lint_path_add(
        self,
        snakefile,
        regex1=re.compile(
            "[a-zA-Z_][a-zA-Z_0-9]* *\\+ *(?P<quote>['\"]).*?/.*?(?P=quote)"
        ),
        regex2=re.compile(
            "(?P<quote>['\"]).*/.*?(?P=quote) *\\+ *{}".format(NAME_PATTERN)
        ),
    ):
        for match in chain(regex1.finditer(snakefile), regex2.finditer(snakefile)):
            line = get_line(match, snakefile)
            yield Lint(
                title="Path composition with '+' in line {}".format(line),
                body='This becomes quickly unreadable. Instead, use pathlib or string formatting with f"...".',
            )

    def lint_envvars(
        self,
        snakefile,
        regex=re.compile("os.environ\[(?P<quote>['\"])(?P<name>.+)?(?P=quote)\]"),
    ):
        for match in regex.finditer(snakefile):
            line = get_line(match, snakefile)
            name = match.group("name")
            if name not in self.workflow.envvars:
                yield Lint(
                    title="Environment variable {} used but not asserted with envvars directive in line {}.".format(
                        name, line
                    ),
                    body="Asserting existence of environment variables with the envvars directive ensures proper error "
                    "messages if the user fails to invoke a workflow with all required environment variables defined. "
                    "Further, it allows snakemake to pass them on in case of distributed execution.",
                    links=[links.envvars],
                )

    def lint_singularity(self, snakefile, regex=re.compile("singularity:")):
        for match in regex.finditer(snakefile):
            line = get_line(match, snakefile)
            yield Lint(
                title="Deprecated singularity directive used for container definition in line {}.".format(
                    line
                ),
                body="Use the container directive instead (it is agnostic of the underlying container runtime).",
                links=[links.containers],
            )


def get_line(match, snakefile):
    return snakefile[: match.start()].count("\n") + 1
