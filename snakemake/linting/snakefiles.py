from collections import nametuple
import re

from snakemake.linting import Linter, Lint

class SnakefileLinter(Linter):
    def item_desc_plain(self, snakefile):
        return "snakefile {}".format(snakefile)
    
    def item_desc_json(self, snakefile):
        return {"snakefile": snakefile}
    
    def lint_absolute_paths(self, snakefile, regex=re.compile("(?P<quote>['\"])/.*?(?P=quote)")):
        for match in regex.finditer(snakefile):
            line = get_line(match, snakefile)
            yield Lint(
                title="Absolute path in line {}.".format(line),
                body="Do not define absolute paths inside of the workflow, since this "
                "renders your workflow irreproducible on other machines. "
                "Use path relative to the working directory instead, or make the path "
                "configurable via a config file.",
                links=["https://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#configuration"]
            )
    
    def lint_mixed_func_and_rules(self, snakefile, rule_regex=re.compile("rule .+?:"), func_regex=re.compile("def .+?:")):
        if rule_regex.search(snakefile) and func_regex.search(snakefile):
            yield Lint(
                title="Mixed rules and functions in same snakefile."
                body="Small one-liner functions used only once should be "
                "defined as lambda expressions. Other functions should be collected "
                "in a common module, e.g. 'rules/common.smk'. This makes the workflow "
                "steps more readable.",
                links=["https://snakemake.readthedocs.io/en/latest/snakefiles/modularization.html#includes"]
            )
    
    def lint_path_add(self, snakefile, regex1=re.compile("[a-zA-Z_][a-zA-Z_0-9]* *+ *(?P<quote>['\"]).*?/.*?(?P=quote)"), regex2=re.compile("(?P<quote>['\"]).*/.*?(?P=quote) *+ *[a-zA-Z_][a-zA-Z_0-9]*")):
        for match in regex.finditer(snakefile):
            line = get_line(match, snakefile)
            yield Lint(
                title="Path composition with '+' in line {}".format(line),
                body="This becomes quickly unreadable. Instead, use pathlib or string formatting with f\"...\"."
            )
    
    def lint_envvars(self, snakefile, regex=re.compile("os.environ\[(?P<quote>['\"])(?P<name>.+)?(?P=quote)\]")):
        for match in regex.finditer(snakefile):
            line = get_line(match)
            name = match.group("name")
            if name not in self.workflow.envvars:
                yield Lint(
                    title="Environment variable {} used but not asserted with envvars directive in line {}.".format(name, line),
                    body="Asserting existence of environment variables with the envvars directive ensures proper error "
                    "messages if the user fails to invoke a workflow with all required environment variables defined. "
                    "Further, it allows snakemake to pass them on in case of distributed execution.",
                    links=["https://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#environment-variables"]
                )

    
def get_line(match, snakefile):
    return snakefile[:match.start()].count("\n") + 1