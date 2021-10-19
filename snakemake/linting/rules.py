from itertools import chain
import re

from snakemake.io import get_wildcard_names, is_flagged
from snakemake.linting import Linter, Lint, links, NAME_PATTERN


class RuleLinter(Linter):
    def item_desc_plain(self, rule):
        return "rule {} (line {}, {})".format(rule.name, rule.lineno, rule.snakefile)

    def item_desc_json(self, rule):
        return {"rule": rule.name, "line": rule.lineno, "snakefile": rule.snakefile}

    def lint_params_prefix(self, rule):
        for param, value in rule.params.items():
            if (
                isinstance(value, str)
                and value
                and any(
                    f.startswith(value)
                    for f in chain(rule.input, rule.output)
                    if isinstance(f, str)
                )
            ):
                yield Lint(
                    title="Param {} is a prefix of input or output file but hardcoded".format(
                        param
                    ),
                    body="If this is meant to represent a file path prefix, it will fail when running "
                    "workflow in environments without a shared filesystem. "
                    "Instead, provide a function that infers the appropriate prefix from the input or "
                    "output file, e.g.: lambda w, input: os.path.splitext(input[0])[0]",
                    links=[links.params, links.input_functions],
                )

    def lint_log_directive(self, rule):
        if not rule.log and not rule.norun and not rule.is_handover:
            yield Lint(
                title="No log directive defined",
                body="Without a log directive, all output will be printed "
                "to the terminal. In distributed environments, this means "
                "that errors are harder to discover. In local environments, "
                "output of concurrent jobs will be mixed and become unreadable.",
                links=[links.log],
            )

    def lint_not_used_params(
        self,
        rule,
        valid_names={
            "input",
            "output",
            "log",
            "params",
            "wildcards",
            "threads",
            "resources",
        },
        regex=re.compile("{{(?P<name>{}).*?}}".format(NAME_PATTERN)),
    ):
        if rule.shellcmd:
            for match in regex.finditer(rule.shellcmd):
                name = match.group("name")

                before = match.start() - 1
                after = match.end()

                if name not in valid_names and (
                    not (before >= 0 and after < len(rule.shellcmd))
                    or (rule.shellcmd[before] != "{" and rule.shellcmd[after] != "}")
                ):
                    yield Lint(
                        title="Shell command directly uses variable {} from outside of the rule".format(
                            name
                        ),
                        body="It is recommended to pass all files as input and output, and non-file parameters "
                        "via the params directive. Otherwise, provenance tracking is less accurate.",
                        links=[links.params],
                    )

    def lint_version(self, rule):
        if rule.version:
            yield Lint(
                title="The version directive is deprecated",
                body="It was meant for documenting tool version, but this has been replaced "
                "by using the conda or container directive.",
                links=[links.package_management, links.containers],
            )

    def lint_dynamic(self, rule):
        for file in chain(rule.output, rule.input):
            if is_flagged(file, "dynamic"):
                yield Lint(
                    title="The dynamic flag is deprecated",
                    body="Use checkpoints instead, which are more powerful and less error-prone.",
                    links=[links.checkpoints],
                )

    def lint_long_run(self, rule):
        func_code = rule.run_func.__code__.co_code

        if rule.is_run and len(func_code) > 70:
            yield Lint(
                title="Migrate long run directives into scripts or notebooks",
                body="Long run directives hamper workflow readability. Use the script or notebook direcive instead. "
                "Note that the script or notebook directive does not involve boilerplate. Similar to run, you "
                "will have direct access to params, input, output, and wildcards."
                "Only use the run direcive for a handful of lines.",
                links=[links.external_scripts, links.notebooks],
            )

    def lint_iofile_by_index(self, rule, regex=re.compile("(input|output)\[[0-9]+\]")):
        if rule.shellcmd and regex.search(rule.shellcmd):
            yield Lint(
                title="Do not access input and output files individually by index in shell commands",
                body="When individual access to input or output files is needed (i.e., just writing '{input}' "
                "is impossible), use names ('{input.somename}') instead of index based access.",
                links=[links.rules],
            )

    def lint_missing_software_definition(self, rule):
        if (
            not rule.norun
            and not rule.is_handover
            and not rule.is_run
            and not (rule.conda_env or rule.container_img)
        ):
            if rule.env_modules:
                yield Lint(
                    title="Additionally specify a conda environment or container for each rule, environment modules are not enough",
                    body="While environment modules allow to document and deploy the required software on a certain "
                    "platform, they lock your workflow in there, disabling easy reproducibility on other machines "
                    "that don't have exactly the same environment modules. Hence env modules (which might be beneficial "
                    "in certain cluster environments), should allways be complemented with equivalent conda "
                    "environments.",
                    links=[links.package_management, links.containers],
                )
            else:
                yield Lint(
                    title="Specify a conda environment or container for each rule.",
                    body="This way, the used software for each specific step is documented, and "
                    "the workflow can be executed on any machine without prerequisites.",
                    links=[links.package_management, links.containers],
                )
