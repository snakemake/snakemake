__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from dataclasses import dataclass
import sys
import textwrap
import tokenize
from typing import Any, Callable, Dict, Generator, List, Optional

import snakemake
from snakemake import common, sourcecache, workflow

dd = textwrap.dedent

INDENT = "\t"


def is_newline(token, newline_tokens=set((tokenize.NEWLINE, tokenize.NL))):
    return token.type in newline_tokens


def is_line_start(token):
    prefix = token.line[: token.start[1]]
    return not prefix or prefix.isspace()


def is_indent(token):
    return token.type == tokenize.INDENT


def is_dedent(token):
    return token.type == tokenize.DEDENT


def is_op(token):
    return token.type == tokenize.OP


def is_greater(token):
    return is_op(token) and token.string == ">"


def is_comma(token):
    return is_op(token) and token.string == ","


def is_name(token):
    return token.type == tokenize.NAME


def is_colon(token):
    return is_op(token) and token.string == ":"


def is_comment(token):
    return token.type == tokenize.COMMENT


def is_string(token):
    return token.type == tokenize.STRING


def is_fstring_start(token):
    return sys.version_info >= (3, 12) and token.type == tokenize.FSTRING_START


def is_eof(token):
    return token.type == tokenize.ENDMARKER


def lineno(token: tokenize.TokenInfo):
    return token.start[0]


class StopAutomaton(Exception):
    def __init__(self, token):
        self.token = token


class TokenAutomaton:
    subautomata: Dict[str, Any] = {}
    deprecated: Dict[str, str] = {}

    def __init__(self, snakefile: "Snakefile", base_indent=0, dedent=0, root=True):
        self.root = root
        self.snakefile = snakefile
        self.state: Callable[[tokenize.TokenInfo], Generator] = None  # type: ignore
        self.base_indent = base_indent
        self.line = 0
        self.indent = 0
        self.was_indented = False
        self.lasttoken = None
        self._dedent = dedent

    @property
    def dedent(self):
        return self._dedent

    @property
    def effective_indent(self):
        return self.base_indent + self.indent - self.dedent

    def indentation(self, token):
        if is_indent(token) or is_dedent(token):
            self.indent = token.end[1] - self.base_indent
            self.was_indented |= self.indent > 0

    def parse_fstring(self, token: tokenize.TokenInfo):
        """
        only for python >= 3.12, since then python changed the
        parsing manner of f-string, see
        [pep-0701](https://peps.python.org/pep-0701)

        Here, we just read where the f-string start and end from tokens.
        Luckily, each token records the content of the line,
        and we can just take what we want there.
        """
        related_lines = token.start[0]
        s = token.line
        isin_fstring = 1
        for t1 in self.snakefile:
            if related_lines < t1.start[0]:
                # go to the next line
                related_lines = t1.start[0]
                s += t1.line
            if t1.type == tokenize.FSTRING_START:
                isin_fstring += 1
            elif t1.type == tokenize.FSTRING_END:
                isin_fstring -= 1
            if isin_fstring == 0:
                break
        # trim those around the f-string
        t = s[token.start[1] : t1.end[1] - len(t1.line)]
        if hasattr(self, "cmd") and self.cmd[-1][1] == token:
            self.cmd[-1] = t, token
        return t

    def consume(self):
        for token in self.snakefile:
            self.indentation(token)
            try:
                for t, orig in self.state(token):
                    # python >= 3.12 only
                    if is_fstring_start(token):
                        t = self.parse_fstring(token)
                    if self.lasttoken == "\n" and not t.isspace():
                        yield INDENT * self.effective_indent, orig
                    yield t, orig
                    self.lasttoken = t
            except tokenize.TokenError as e:
                self.error(
                    str(e).split(",")[0].strip("()''"), token
                )  # TODO the inferred line number seems to be wrong sometimes

    def error(self, msg, token, naming_hint=None):
        if naming_hint is not None:
            msg += (
                f" The keyword {naming_hint} has a special meaning in Snakemake. "
                "If you named a variable or function like this, please rename it to "
                "avoid the conflict."
            )
        raise SyntaxError(msg, (self.snakefile.path, lineno(token), None, None))

    def subautomaton(self, automaton, *args, token=None, **kwargs):
        if automaton in self.deprecated:
            assert (
                token is not None
            ), "bug: deprecation encountered but subautomaton not called with a token"
            self.error(
                f"Keyword {automaton} is deprecated. {self.deprecated[automaton]}",
                token,
            )
        return self.subautomata[automaton](
            self.snakefile,
            *args,
            base_indent=self.base_indent + self.indent,
            dedent=self.dedent,
            root=False,
            **kwargs,
        )


class KeywordState(TokenAutomaton):
    prefix = ""
    start: Callable[[], Generator[str, None, None]]

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile, base_indent=base_indent, dedent=dedent, root=root)
        self.line = 0
        self.state = self.colon

    @property
    def keyword(self):
        return self.__class__.__name__.lower()[len(self.prefix) :]

    def end(self):
        # Add newline to prevent https://github.com/snakemake/snakemake/issues/1943
        yield "\n"
        yield ")"

    def decorate_end(self, token):
        for t in self.end():
            if isinstance(t, tuple):
                yield t
            else:
                yield t, token

    def colon(self, token):
        if is_colon(token):
            self.state = self.block
            for t in self.start():
                yield t, token
        else:
            self.error(f"Colon expected after keyword {self.keyword}.", token)

    def is_block_end(self, token):
        return (self.line and self.indent <= 0) or is_eof(token)

    def block(self, token, force_block_end=False):
        if self.lasttoken == "\n" and is_comment(token):
            # ignore lines containing only comments
            self.line -= 1
        if force_block_end or self.is_block_end(token):
            yield from self.decorate_end(token)
            yield "\n", token
            raise StopAutomaton(token)

        if is_newline(token):
            self.line += 1
            yield token.string, token
        elif not (is_indent(token) or is_dedent(token)):
            if is_comment(token):
                yield token.string, token
            else:
                yield from self.block_content(token)

    def yield_indent(self, token):
        return token.string, token

    def block_content(self, token):
        yield token.string, token


class GlobalKeywordState(KeywordState):
    def start(self):
        yield f"workflow.{self.keyword}("


class DecoratorKeywordState(KeywordState):
    decorator: Optional[str] = None
    args: List[str] = []

    def start(self):
        yield f"@workflow.{self.decorator}"
        yield "\n"
        yield "def __{}({}):".format(self.decorator, ", ".join(self.args))

    def end(self):
        yield ""


class RuleKeywordState(KeywordState):
    def __init__(self, snakefile, base_indent=0, dedent=0, root=True, rulename=None):
        super().__init__(snakefile, base_indent=base_indent, dedent=dedent, root=root)
        self.rulename = rulename

    def start(self):
        yield "\n"
        yield f"@workflow.{self.keyword}("


class SectionKeywordState(KeywordState):
    def start(self):
        yield f", {self.keyword}="

    def end(self):
        # no end needed
        return list()


# Global keyword states


class Envvars(GlobalKeywordState):
    @property
    def keyword(self):
        return "register_envvars"


class Include(GlobalKeywordState):
    pass


class Workdir(GlobalKeywordState):
    pass


class Configfile(GlobalKeywordState):
    pass


# PEPs


class Pepfile(GlobalKeywordState):
    @property
    def keyword(self):
        return "set_pepfile"


class Pepschema(GlobalKeywordState):
    pass


class Report(GlobalKeywordState):
    pass


class Scattergather(GlobalKeywordState):
    pass


class Storage(GlobalKeywordState):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tag = None
        self.state = self.register_tag

    def start(self):
        yield f"workflow.storage_registry.register_storage(tag={self.tag!r}, "

    def register_tag(self, token):
        if is_name(token):
            self.tag = token.string
        elif is_colon(token):
            self.state = self.block
            for t in self.start():
                yield t, token
        else:
            self.error(
                "Expected name or colon after storage keyword.",
                token,
                naming_hint="storage",
            )


class ResourceScope(GlobalKeywordState):
    err_msg = (
        "Invalid scope: {resource}={scope}. Scope must be set to either 'local' or "
        "'global'"
    )
    current_resource = ""

    def block_content(self, token):
        if is_name(token):
            self.current_resource = token.string
        if is_string(token) and self.lasttoken == "=":
            if token.string[1:][:-1] not in ["local", "global"]:
                self.error(
                    self.err_msg.format(
                        resource=self.current_resource, scope=token.string
                    ),
                    token,
                )
        yield token.string, token


class Ruleorder(GlobalKeywordState):
    def block_content(self, token):
        if is_greater(token):
            yield ",", token
        elif is_name(token):
            yield repr(token.string), token
        else:
            self.error(
                "Expected a descending order of rule names, "
                "e.g. rule1 > rule2 > rule3 ...",
                token,
            )


class GlobalWildcardConstraints(GlobalKeywordState):
    @property
    def keyword(self):
        return "global_wildcard_constraints"


class GlobalSingularity(GlobalKeywordState):
    @property
    def keyword(self):
        return "global_container"


class GlobalContainer(GlobalKeywordState):
    @property
    def keyword(self):
        return "global_container"


class GlobalContainerized(GlobalKeywordState):
    @property
    def keyword(self):
        return "global_containerized"


class GlobalConda(GlobalKeywordState):
    @property
    def keyword(self):
        return "global_conda"


class Localrules(GlobalKeywordState):
    def block_content(self, token):
        if is_comma(token):
            yield ",", token
        elif is_name(token):
            yield repr(token.string), token
        else:
            self.error(
                "Expected a comma separated list of rules that shall "
                "not be executed by the cluster command.",
                token,
            )


# Rule keyword states


class Name(RuleKeywordState):
    pass


class Input(RuleKeywordState):
    pass


class Output(RuleKeywordState):
    pass


class Params(RuleKeywordState):
    pass


class Threads(RuleKeywordState):
    pass


class Retries(RuleKeywordState):
    pass


class Shadow(RuleKeywordState):
    pass


class Resources(RuleKeywordState):
    pass


class Priority(RuleKeywordState):
    pass


class Version(RuleKeywordState):
    pass


class Log(RuleKeywordState):
    pass


class Message(RuleKeywordState):
    pass


class Benchmark(RuleKeywordState):
    pass


class Conda(RuleKeywordState):
    pass


class Singularity(RuleKeywordState):
    @property
    def keyword(self):
        return "container"


class Container(RuleKeywordState):
    pass


class Containerized(RuleKeywordState):
    pass


class EnvModules(RuleKeywordState):
    pass


class Group(RuleKeywordState):
    pass


class Cache(RuleKeywordState):
    @property
    def keyword(self):
        return "cache_rule"


class DefaultTarget(RuleKeywordState):
    @property
    def keyword(self):
        return "default_target_rule"


class Handover(RuleKeywordState):
    pass


class WildcardConstraints(RuleKeywordState):
    @property
    def keyword(self):
        return "register_wildcard_constraints"


class LocalRule(RuleKeywordState):
    pass


class Run(RuleKeywordState):
    def __init__(self, snakefile, rulename, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile, base_indent=base_indent, dedent=dedent, root=root)
        self.rulename = rulename
        self.content = 0

    def start(self):
        yield "@workflow.run"
        yield "\n"
        yield (
            "def __rule_{rulename}(input, output, params, wildcards, threads, "
            "resources, log, rule, conda_env, container_img, "
            "singularity_args, use_singularity, env_modules, bench_record, jobid, "
            "is_shell, bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, "
            "conda_base_path, basedir, sourcecache_path, runtime_sourcecache_path, {rule_func_marker}=True):".format(
                rulename=(
                    self.rulename
                    if self.rulename is not None
                    else self.snakefile.rulecount
                ),
                rule_func_marker=common.RULEFUNC_CONTEXT_MARKER,
            )
        )

    def end(self):
        yield ""

    def block_content(self, token):
        self.content += 1
        yield token.string, token

    def is_block_end(self, token):
        return (self.content and self.line and self.indent <= 0) or is_eof(token)


class AbstractCmd(Run):
    overwrite_cmd: Optional[str] = None
    start_func: Optional[str] = None
    end_func: Optional[str] = None

    def __init__(self, snakefile, rulename, base_indent=0, dedent=0, root=True):
        super().__init__(
            snakefile, rulename, base_indent=base_indent, dedent=dedent, root=root
        )
        self.cmd: list[tuple[str, tokenize.TokenInfo]] = []
        self.token = None
        if self.overwrite_cmd is not None:
            self.block_content = self.overwrite_block_content  # type: ignore

    def is_block_end(self, token):
        return (self.line and self.indent <= 0) or is_eof(token)

    def start(self):
        if self.start_func is not None:
            yield self.start_func
            yield "("

    def args(self):
        yield from []

    def end(self):
        # the end is detected. So we can safely reset the indent to zero here
        self.indent = 0
        yield "\n"
        yield ")"
        yield "\n"
        for t in super().start():
            yield t
        yield "\n"
        yield INDENT * (self.effective_indent + 1)
        yield self.end_func
        yield "("
        yield from self.cmd
        yield from self.args()
        yield "\n"
        yield ")"
        for t in super().end():
            yield t

    def decorate_end(self, token):
        if self.token is None:
            # no block after shell keyword
            self.error(
                "Command must be given as string after the shell keyword.", token
            )
        yield from super().decorate_end(self.token)

    def block_content(self, token):
        self.token = token
        self.cmd.append((token.string, token))
        yield token.string, token

    def overwrite_block_content(self, token):
        if self.token is None:
            self.token = token
            cmd = repr(self.overwrite_cmd)
            self.cmd.append((cmd, token))
            yield cmd, token


class Shell(AbstractCmd):
    start_func = "@workflow.shellcmd"
    end_func = "shell"

    def args(self):
        yield ", bench_record=bench_record, bench_iteration=bench_iteration"


class Script(AbstractCmd):
    start_func = "@workflow.script"
    end_func = "script"

    def args(self):
        yield (
            ", basedir, input, output, params, wildcards, threads, resources, log, "
            "config, rule, conda_env, conda_base_path, container_img, singularity_args, env_modules, "
            "bench_record, jobid, bench_iteration, cleanup_scripts, shadow_dir, sourcecache_path, "
            "runtime_sourcecache_path"
        )


class Notebook(Script):
    start_func = "@workflow.notebook"
    end_func = "notebook"

    def args(self):
        yield (
            ", basedir, input, output, params, wildcards, threads, resources, log, "
            "config, rule, conda_env, conda_base_path, container_img, singularity_args, env_modules, "
            "bench_record, jobid, bench_iteration, cleanup_scripts, shadow_dir, "
            "edit_notebook, sourcecache_path, runtime_sourcecache_path"
        )


class Wrapper(Script):
    start_func = "@workflow.wrapper"
    end_func = "wrapper"

    def args(self):
        yield (
            ", input, output, params, wildcards, threads, resources, log, "
            "config, rule, conda_env, conda_base_path, container_img, singularity_args, env_modules, "
            "bench_record, workflow.workflow_settings.wrapper_prefix, jobid, bench_iteration, "
            "cleanup_scripts, shadow_dir, sourcecache_path, runtime_sourcecache_path"
        )


class TemplateEngine(Script):
    start_func = "@workflow.template_engine"
    end_func = "render_template"

    def args(self):
        yield (", input, output, params, wildcards, config, rule")


class CWL(Script):
    start_func = "@workflow.cwl"
    end_func = "cwl"

    def args(self):
        yield (
            ", basedir, input, output, params, wildcards, threads, resources, log, "
            "config, rule, use_singularity, bench_record, jobid, sourcecache_path, "
            "runtime_sourcecache_path"
        )


rule_property_subautomata = dict(
    name=Name,
    input=Input,
    output=Output,
    params=Params,
    threads=Threads,
    resources=Resources,
    retries=Retries,
    priority=Priority,
    log=Log,
    message=Message,
    benchmark=Benchmark,
    conda=Conda,
    singularity=Singularity,
    container=Container,
    containerized=Containerized,
    envmodules=EnvModules,
    wildcard_constraints=WildcardConstraints,
    shadow=Shadow,
    group=Group,
    cache=Cache,
    handover=Handover,
    default_target=DefaultTarget,
    localrule=LocalRule,
)
rule_property_deprecated = dict(
    version="Use conda or container directive instead (see docs)."
)


class Rule(GlobalKeywordState):
    subautomata = dict(
        run=Run,
        shell=Shell,
        script=Script,
        notebook=Notebook,
        wrapper=Wrapper,
        template_engine=TemplateEngine,
        cwl=CWL,
        **rule_property_subautomata,
    )
    deprecated = rule_property_deprecated

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile, base_indent=base_indent, dedent=dedent, root=root)
        self.state = self.name
        self.lineno = None
        self.rulename = None
        self.run = False
        self.snakefile.rulecount += 1

    def start(self, aux=""):
        yield (
            f"@workflow.rule(name={self.rulename!r}, lineno={self.lineno}, "
            f"snakefile={self.snakefile.path!r}{aux})"
        )

    def end(self):
        if not self.run:
            yield "@workflow.norun()"
            yield "\n"
            for t in self.subautomaton("run", rulename=self.rulename).start():
                yield t
            # the end is detected.
            # So we can safely reset the indent to zero here
            self.indent = 0
            yield "\n"
            yield INDENT * (self.effective_indent + 1)
            yield "pass"

    def name(self, token):
        if is_name(token):
            self.rulename = token.string
        elif is_colon(token):
            self.lineno = self.snakefile.lines + 1
            self.state = self.block
            for t in self.start():
                yield t, token
        else:
            self.error(
                "Expected name or colon after rule or checkpoint keyword.",
                token,
                naming_hint="rule",
            )

    def block_content(self, token):
        if is_name(token):
            try:
                if (
                    token.string == "run"
                    or token.string == "shell"
                    or token.string == "script"
                    or token.string == "wrapper"
                    or token.string == "notebook"
                    or token.string == "template_engine"
                    or token.string == "cwl"
                ):
                    if self.run:
                        raise self.error(
                            "Multiple run/shell/script/notebook/wrapper/template_engine/cwl "
                            "keywords in rule {}.".format(self.rulename),
                            token,
                        )
                    self.run = True
                elif self.run:
                    raise self.error(
                        "No rule keywords allowed after "
                        "run/shell/script/notebook/wrapper/template_engine/cwl in "
                        "rule {}.".format(self.rulename),
                        token,
                    )
                for t in self.subautomaton(
                    token.string, token=token, rulename=self.rulename
                ).consume():
                    yield t
            except KeyError:
                self.error(
                    f"Unexpected keyword {token.string} in rule definition",
                    token,
                )
            except StopAutomaton as e:
                self.indentation(e.token)
                for t in self.block(e.token):
                    yield t
        elif is_comment(token):
            yield "\n", token
            yield token.string, token
        elif is_string(token):
            yield "\n", token
            yield f"@workflow.docstring({token.string})", token
        else:
            self.error(
                "Expecting rule keyword, comment or docstrings "
                "inside a rule definition.",
                token,
            )

    @property
    def dedent(self):
        return self.indent


class Checkpoint(Rule):
    def start(self):
        yield from super().start(aux=", checkpoint=True")


class OnSuccess(DecoratorKeywordState):
    decorator = "onsuccess"
    args = ["log"]


class OnError(DecoratorKeywordState):
    decorator = "onerror"
    args = ["log"]


class OnStart(DecoratorKeywordState):
    decorator = "onstart"
    args = ["log"]


# modules


class ModuleKeywordState(SectionKeywordState):
    prefix = "Module"


class ModuleSnakefile(ModuleKeywordState):
    pass


class ModulePrefix(ModuleKeywordState):
    pass


class ModuleMetaWrapper(ModuleKeywordState):
    @property
    def keyword(self):
        return "meta_wrapper"


class ModuleConfig(ModuleKeywordState):
    pass


class ModuleSkipValidation(ModuleKeywordState):
    @property
    def keyword(self):
        return "skip_validation"


class ModuleReplacePrefix(ModuleKeywordState):
    @property
    def keyword(self):
        return "replace_prefix"


class Module(GlobalKeywordState):
    subautomata = dict(
        snakefile=ModuleSnakefile,
        meta_wrapper=ModuleMetaWrapper,
        config=ModuleConfig,
        skip_validation=ModuleSkipValidation,
        replace_prefix=ModuleReplacePrefix,
        prefix=ModulePrefix,
    )

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile, base_indent=base_indent, dedent=dedent, root=root)
        self.state = self.name
        self.has_snakefile = False
        self.has_meta_wrapper = False
        self.has_name = False
        self.primary_token = None

    def end(self):
        if not (self.has_snakefile or self.has_meta_wrapper):
            self.error(
                "A module needs either a path to a Snakefile or a meta wrapper URL.",
                self.primary_token,
            )
        yield ")"

    def name(self, token):
        if is_name(token):
            yield f"workflow.module({token.string!r}", token
            self.has_name = True
        elif is_colon(token) and self.has_name:
            self.primary_token = token
            self.state = self.block
        else:
            self.error(
                "Expected name after module keyword.", token, naming_hint="module"
            )

    def block_content(self, token):
        if is_name(token):
            try:
                if token.string == "snakefile":
                    self.has_snakefile = True
                if token.string == "meta_wrapper":
                    self.has_meta_wrapper = True
                for t in self.subautomaton(token.string, token=token).consume():
                    yield t
            except KeyError:
                self.error(
                    "Unexpected keyword {} in "
                    "module definition".format(token.string),
                    token,
                )
            except StopAutomaton as e:
                self.indentation(e.token)
                for t in self.block(e.token):
                    yield t
        elif is_comment(token):
            yield "\n", token
            yield token.string, token
        elif is_string(token):
            # ignore docstring
            pass
        else:
            self.error(
                "Expecting module keyword, comment or docstrings "
                "inside a module definition.",
                token,
            )


class UseRule(GlobalKeywordState):
    subautomata = rule_property_subautomata
    deprecated = rule_property_deprecated

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile, base_indent=base_indent, dedent=dedent, root=root)
        self.state = self.state_keyword_rule
        self.rules = []
        self.exclude_rules = []
        self.has_with = False
        self.name_modifier = []
        self.from_module = None
        self._with_block = []
        self.lineno = self.snakefile.lines + 1

    def end(self):
        name_modifier = "".join(self.name_modifier) if self.name_modifier else None
        yield "@workflow.userule(rules={!r}, from_module={!r}, exclude_rules={!r}, name_modifier={!r}, lineno={})".format(
            self.rules, self.from_module, self.exclude_rules, name_modifier, self.lineno
        )
        yield "\n"

        if self._with_block:
            # yield with block
            yield from self._with_block

            yield "@workflow.run"
            yield "\n"

        rulename = self.rules[0]
        if rulename == "*":
            rulename = "__allrules__"
        yield f"def __userule_{self.from_module}_{rulename}():"
        # the end is detected.
        # So we can safely reset the indent to zero here
        self.indent = 0
        yield "\n"
        yield INDENT * (self.effective_indent + 1)
        yield "pass"

    def state_keyword_rule(self, token):
        if is_name(token) and token.string == "rule":
            self.state = self.state_rules_rule
            yield from ()
        else:
            self.error("Expecting keyword 'rule' after keyword 'use'", token)

    def state_rules_rule(self, token):
        if is_name(token):
            if token.string == "from" or token.string == "as" and not self.rules:
                self.error("Expecting rule names after 'use rule' statement.", token)

            self.rules.append(token.string)
            self.state = self.state_rules_comma_or_end
            yield from ()
        elif is_op(token):
            if token.string == "*":
                self.rules.append(token.string)
                self.state = self.state_rules_end
                yield from ()
            else:
                self.error(
                    "Expecting rule name or '*' after 'use rule' statement.", token
                )
        else:
            self.error(
                "Expecting rule listing (comma separated) after 'use rule' statement.",
                token,
            )
        # TODO newline and parentheses handling

    def state_rules_end(self, token):
        if is_name(token) and token.string == "from":
            self.state = self.state_from
            yield from ()
        else:
            self.error(
                "Expecting list of rules in 'use rule' statement to end with keyword 'from'.",
                token,
            )

    def state_rules_comma_or_end(self, token):
        if is_name(token):
            if token.string == "from" or token.string == "as":
                if not self.rules:
                    self.error(
                        "Expecting rule names after 'use rule' statement.", token
                    )
                if token.string == "from":
                    self.state = self.state_from
                else:
                    self.state = self.state_as
                yield from ()
            else:
                self.error(
                    "Expecting list of rules in 'use rule' statement to end with keyword 'from'.",
                    token,
                )
        elif is_comma(token):
            self.state = self.state_rules_rule
            yield from ()
        else:
            self.error(
                "Unexpected token in list of rules within 'use rule' statement.", token
            )

    def state_from(self, token):
        if is_name(token):
            self.state = self.state_modifier
            self.from_module = token.string
            yield from ()
        else:
            self.error(
                "Expecting module name after 'from' keyword in 'use rule' statement.",
                token,
            )

    def state_modifier(self, token):
        if is_name(token):
            if token.string == "as" and not self.name_modifier:
                self.state = self.state_as
                yield from ()
            elif token.string == "exclude":
                self.state = self.state_exclude
                yield from ()
            elif token.string == "with":
                yield from self.handle_with(token)
            else:
                self.error(
                    "Expecting at most one 'as' or 'with' statement, or the end of the line.",
                    token,
                )
        elif is_newline(token) or is_comment(token) or is_eof(token):
            # end of the statement, close block manually
            yield from self.block(token, force_block_end=True)
        else:
            self.error(
                "Expecting either 'as', 'with' or end of line in 'use rule' statement.",
                token,
            )

    def handle_with(self, token):
        if "*" in self.rules:
            self.error(
                "Keyword 'with' in 'use rule' statement is not allowed in combination with rule pattern '*'.",
                token,
            )
        self.has_with = True
        self.state = self.state_with
        yield from ()

    def state_as(self, token):
        if is_name(token):
            if token.string != "with":
                self.name_modifier.append(token.string)
                yield from ()
            else:
                yield from self.handle_with(token)
        elif is_op(token) and token.string == "*":
            self.name_modifier.append(token.string)
            yield from ()
        elif is_newline(token) or is_comment(token) or is_eof(token):
            # end of the statement, close block manually
            yield from self.block(token, force_block_end=True)
        else:
            self.error(
                "Expecting rulename modifying pattern (e.g. modulename_*) after 'as' keyword.",
                token,
            )

    def state_with(self, token):
        if is_colon(token):
            self.state = self.block
            yield from ()
        else:
            self.error(
                "Expecting colon after 'with' keyword in 'use rule' statement.", token
            )

    def state_exclude(self, token):
        if is_name(token):
            self.exclude_rules.append(token.string)
            self.state = self.state_exclude_comma_or_end
            yield from ()
        else:
            self.error(
                "Expecting rule name(s) after 'exclude' keyword in 'use rule' statement.",
                token,
            )

    def state_exclude_comma_or_end(self, token):
        if is_name(token):
            if token.string == "from" or token.string == "as":
                if not self.exclude_rules:
                    self.error("Expecting rule names after 'exclude' statement.", token)
                if token.string == "from":
                    self.state = self.state_from
                else:
                    self.state = self.state_as
                yield from ()
            else:
                yield from ()
        elif is_comma(token):
            self.state = self.state_exclude
            yield from ()
        else:
            self.state = self.state_modifier
            yield from ()

    def block_content(self, token):
        if is_comment(token):
            yield "\n", token
            yield token.string, token
        elif is_name(token):
            try:
                self._with_block.extend(
                    self.subautomaton(token.string, token=token).consume()
                )
                yield from ()
            except KeyError:
                self.error(
                    f"Unexpected keyword {token.string} in rule definition",
                    token,
                )
            except StopAutomaton as e:
                self.indentation(e.token)
                yield from self.block(e.token)
        else:
            self.error(
                "Expecting a keyword or comment "
                "inside a 'use rule ... with:' statement.",
                token,
            )

    @property
    def dedent(self):
        return self.indent


class Python(TokenAutomaton):
    subautomata = dict(
        envvars=Envvars,
        include=Include,
        workdir=Workdir,
        configfile=Configfile,
        pepfile=Pepfile,
        pepschema=Pepschema,
        report=Report,
        ruleorder=Ruleorder,
        rule=Rule,
        checkpoint=Checkpoint,
        localrules=Localrules,
        onsuccess=OnSuccess,
        onerror=OnError,
        onstart=OnStart,
        wildcard_constraints=GlobalWildcardConstraints,
        singularity=GlobalSingularity,
        container=GlobalContainer,
        containerized=GlobalContainerized,
        conda=GlobalConda,
        scattergather=Scattergather,
        storage=Storage,
        resource_scopes=ResourceScope,
        module=Module,
        use=UseRule,
    )
    deprecated = dict(subworkflow="Use module directive instead (see docs).")

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile, base_indent=base_indent, dedent=dedent, root=root)
        self.state = self.python

    def python(self, token: tokenize.TokenInfo):
        if not (is_indent(token) or is_dedent(token)):
            if (
                self.lasttoken is None
                or self.lasttoken.isspace()
                and is_line_start(token)
            ):
                try:
                    for t in self.subautomaton(token.string, token=token).consume():
                        yield t
                except KeyError:
                    yield token.string, token
                except StopAutomaton as e:
                    self.indentation(e.token)
                    for t in self.python(e.token):
                        yield t
            else:
                yield token.string, token


class Snakefile:
    def __init__(
        self,
        path: "sourcecache.SourceFile",
        workflow: "workflow.Workflow",
        rulecount=0,
    ):
        self.path = path.get_path_or_uri()
        self.file = workflow.sourcecache.open(path)
        self.tokens = tokenize.generate_tokens(self.file.readline)
        self.rulecount = rulecount
        self.lines = 0

    def __next__(self):
        return next(self.tokens)

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


def format_tokens(tokens) -> Generator[str, None, None]:
    t_: Optional[str] = None
    for t in tokens:
        if t_ and not t.isspace() and not t_.isspace():
            yield " "
        yield t
        t_ = t


def parse(
    path,
    workflow: "workflow.Workflow",
    linemap: Dict[int, int],
    overwrite_shellcmd=None,
    rulecount=0,
):
    Shell.overwrite_cmd = overwrite_shellcmd
    with Snakefile(path, workflow, rulecount=rulecount) as snakefile:
        automaton = Python(snakefile)
        compilation = list()
        for t, orig_token in automaton.consume():
            line_number = lineno(orig_token)
            linemap |= {
                i: line_number
                for i in range(snakefile.lines + 1, snakefile.lines + t.count("\n") + 1)
            }
            snakefile.lines += t.count("\n")
            compilation.append(t)
    join_compilation = "".join(format_tokens(compilation))
    if linemap:
        last = max(linemap)
        linemap[last + 1] = linemap[last]
    return join_compilation, snakefile.rulecount
