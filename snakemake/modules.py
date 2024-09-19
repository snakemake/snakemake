__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
import types
import re
from copy import copy

from .exceptions import WorkflowError
from .path_modifier import PathModifier
from . import wrapper, ruleinfo
from . import workflow as _workflow
from . import rules as _rules


def _default_modify_rulename(rulename: str):
    return rulename


def get_name_modifier_func(
    rules: list[str] | None = None,
    name_modifier: str | None = None,
    parent_modifier: "WorkflowModifier | None" = None,
):
    if name_modifier is None:
        return _default_modify_rulename
    if parent_modifier is None:
        parent_modifier_func = _default_modify_rulename
    else:
        parent_modifier_func = parent_modifier.modify_rulename
    if "*" in name_modifier:
        return lambda rulename: parent_modifier_func(
            name_modifier.replace("*", rulename)
        )
    if rules and len(rules) > 1:
        raise SyntaxError(
            "Multiple rules in 'use rule' statement but name modification ('as' statement) does not contain a wildcard '*'."
        )
    return lambda rulename: parent_modifier_func(name_modifier)


class ModuleInfo:
    def __init__(
        self,
        workflow: "_workflow.Workflow",
        name,
        snakefile: str | None = None,
        meta_wrapper=None,
        config: dict | None = None,
        skip_validation=False,
        replace_prefix=None,
        prefix=None,
    ):
        self.workflow = workflow
        self.name = name
        self.snakefile = snakefile
        self.meta_wrapper = meta_wrapper
        self.config = config
        self.skip_validation = skip_validation
        self.parent_modifier = self.workflow.modifier
        self.rule_proxies = _rules.Rules()
        self.namespace = types.ModuleType(name)
        self.stack_len = 0

        if prefix is not None:
            if isinstance(prefix, Path):
                prefix = str(prefix)
            if not isinstance(prefix, str):
                raise WorkflowError(
                    "Prefix definition in module statement must be string or Path."
                )
            if replace_prefix is not None:
                raise WorkflowError(
                    "Module definition contains both prefix and replace_prefix. "
                    "Only one at a time is allowed."
                )

        self.replace_prefix = replace_prefix
        self.prefix = prefix

    def load_ruleinfos(self, skip_global_report_caption=False):
        exclude_rules = None
        name_modifier = None
        ruleinfo = None
        rules = ["*"]
        snakefile = self.get_snakefile()
        self.modifier = WorkflowModifier(
            self.workflow,
            config=self.config,
            base_snakefile=snakefile,
            skip_configfile=self.config is not None,
            skip_validation=self.skip_validation,
            skip_global_report_caption=skip_global_report_caption,
            rule_exclude_list=exclude_rules,
            rule_whitelist=[],
            resolved_rulename_modifier=get_name_modifier_func(
                rules, name_modifier, parent_modifier=self.parent_modifier
            ),
            local_rulename_modifier=get_name_modifier_func(rules, name_modifier),
            ruleinfo_overwrite=ruleinfo,
            allow_rule_overwrite=True,
            namespace=self.name,
            replace_prefix=self.replace_prefix,
            prefix=self.prefix,
            replace_wrapper_tag=self.get_wrapper_tag(),
            rule_proxies=self.rule_proxies,
        )
        workflow_rules, self.workflow._rules = self.workflow._rules, {}
        with self.modifier:
            self.workflow.include(snakefile, overwrite_default_target=True)
            self.stack_len = len(self.workflow.included_stack) + 1
        self.workflow._rules = workflow_rules
        if self.name:
            # update globals in inner modules
            self.namespace.__dict__.update(self.modifier.globals)
            self.workflow.globals[self.modifier.namespace] = self.namespace
            for i in self.rule_proxies._cache_rules.values():
                i[0].rule_whitelist = None

    def use_rules(
        self,
        rules: list[str],
        name_modifier=None,
        exclude_rules=None,
        ruleinfo: "_workflow.RuleInfo | None" = None,
    ):
        """
        Rules are selected from `self.rule_proxies._cache_rules`,
            the `ruleinfo_overwrite` will be update,
            and adjusted to the correct order as `Workflow.include`.
        Finally,
            avail rules and ruleorder will be parsed to the parent_modifier.
        """
        name_modifier_func = get_name_modifier_func(
            None, name_modifier, parent_modifier=self.parent_modifier
        )
        # ideally, _old_ruleinfo is None
        _old_ruleinfo, self.modifier.ruleinfo_overwrite = (
            self.modifier.ruleinfo_overwrite,
            ruleinfo,
        )
        snakefile, stacks = self.avail_rules_stacks(
            self.get_rule_whitelist(rules), exclude_rules
        )
        self._include(snakefile, stacks, name_modifier_func, True)
        with self.modifier:
            self.parent_modifier.inherit_rule_proxies(self.modifier, self.stack_len)
            self.parent_modifier.inherit_ruleorder(self.modifier, name_modifier_func)
        self.modifier.ruleinfo_overwrite = _old_ruleinfo

    def get_snakefile(self):
        if self.meta_wrapper:
            return wrapper.get_path(
                self.meta_wrapper + "/test/Snakefile",
                self.workflow.workflow_settings.wrapper_prefix,
            )
        elif self.snakefile:
            return self.snakefile
        else:
            raise WorkflowError(
                "Module statement must either define snakefile or meta_wrapper to use."
            )

    def get_wrapper_tag(self):
        if self.meta_wrapper:
            if wrapper.is_url(self.meta_wrapper):
                raise WorkflowError(
                    "meta_wrapper directive of module statement currently does not support full URLs."
                )
            return self.meta_wrapper.split("/", 1)[0]
        return None

    def get_rule_whitelist(self, rules: list[str]):
        if "*" in rules:
            if len(rules) != 1:
                raise SyntaxError(
                    "The 'use rule' statement uses a wildcard '*' but lists multiple additional rules."
                )
            else:
                return None
        return set(rules)

    def avail_rules_stacks(
        self,
        rule_whitelist: set[str] | None = None,
        exclude_rules: list[str] | None = None,
    ):
        """
        Here, rules are reordered back as loaded from the module.
        It is used to inherent correct default rule from the module.
        """
        excludes = set() if exclude_rules is None else set(exclude_rules)
        rulenames = (
            self.rule_proxies._cache_rules if rule_whitelist is None else rule_whitelist
        )
        pseudo_stacks: list[tuple[set[str | None], list[str | tuple]]] = [(set(), [])]
        last_stack_key: tuple[str | None, int] = None, self.stack_len
        for rulename in rulenames:
            if rulename in excludes:
                continue
            modifier, ruleinfo_, lineno, snakefile, checkpoint, stack_len = (
                self.rule_proxies._cache_rules[rulename]
            )
            if last_stack_key != (snakefile, stack_len):
                if not (
                    last_stack_key == (None, self.stack_len)
                    and stack_len == self.stack_len
                ):
                    if last_stack_key[1] > stack_len:
                        # step out from last include
                        for _p in range(stack_len, last_stack_key[1]):
                            pseudo_stacks[-2][1].append(pseudo_stacks.pop(-1))
                    elif last_stack_key[1] == stack_len:
                        # a new include
                        pseudo_stacks[-2][1].append(pseudo_stacks.pop(-1))
                        pseudo_stacks.append((set(), []))
                    elif last_stack_key[1] < stack_len:
                        for _p in range(last_stack_key[1], stack_len):
                            pseudo_stacks.append((set(), []))
                    else:
                        raise NotImplementedError
                last_stack_key = snakefile, stack_len
            pseudo_stacks[-1][0].add(snakefile)
            pseudo_stacks[-1][1].append(rulename)
        while len(pseudo_stacks) > 1:
            pseudo_stacks[-2][1].append(pseudo_stacks.pop(-1))
        if len(pseudo_stacks[0][0]) > 0:
            (included,) = pseudo_stacks[0][0]
        else:
            # keep consistent only, we know it isn't
            included = self.snakefile
        return included, pseudo_stacks[0][1]

    def _include(
        self,
        snakefile: str | None,
        rule_stacks: list[str | tuple],
        name_modifier_func=_default_modify_rulename,
        overwrite_default_target=False,
    ):
        with _workflow._IncludeWrapper(
            self.workflow,
            snakefile,
            overwrite_default_target=overwrite_default_target,
            module_use=True,
        ):
            for rule_s in rule_stacks:
                if isinstance(rule_s, tuple):
                    if len(rule_s[0]) > 0:
                        (included,) = rule_s[0]
                    else:
                        # keep consistent only, we know it isn't
                        included = snakefile
                    self._include(included, rule_s[1], name_modifier_func)
                else:
                    resolved_rulename = name_modifier_func(rule_s)
                    modifier, ruleinfo_, lineno, snakefile, checkpoint, stack_len = (
                        self.rule_proxies._cache_rules[rule_s]
                    )
                    orig_ruleinfo = copy(ruleinfo_)
                    with modifier:
                        self.workflow.rule(
                            name=resolved_rulename,
                            lineno=lineno,
                            snakefile=snakefile,
                            checkpoint=checkpoint,
                        )(orig_ruleinfo)


class WorkflowModifier:

    def __init__(
        self,
        workflow: "_workflow.Workflow",
        parent_modifier: "WorkflowModifier | None" = None,
        globals=None,
        config=None,
        base_snakefile=None,
        skip_configfile=False,
        skip_validation=False,
        skip_global_report_caption=False,
        resolved_rulename_modifier=_default_modify_rulename,
        local_rulename_modifier=_default_modify_rulename,
        rule_whitelist=None,
        rule_exclude_list=None,
        ruleinfo_overwrite: "ruleinfo.RuleInfo | None" = None,
        allow_rule_overwrite=False,
        replace_prefix=None,
        prefix=None,
        replace_wrapper_tag=None,
        namespace=None,
        rule_proxies: "_rules.Rules | None" = None,
    ):
        if parent_modifier is None:
            # default settings for globals if not inheriting from parent
            self.globals = (
                globals if globals is not None else dict(workflow.vanilla_globals)
            )
            self.wildcard_constraints: dict[str, str] = dict()
            self.rules: set["_rules.Rule"] = set()
            self.globals["rules"] = rule_proxies or _rules.Rules()
            if config is not None:
                self.globals["config"] = config
        else:
            # init with values from parent modifier
            self.globals = parent_modifier.globals
            self.wildcard_constraints = parent_modifier.wildcard_constraints
            self.rules = parent_modifier.rules
            if config is not None:
                raise WorkflowError(
                    "module cannot be loaded with different config with parent modifier exists.\n"
                )

        self.workflow = workflow
        self.base_snakefile = base_snakefile

        self.skip_configfile = skip_configfile
        self.resolved_rulename_modifier = resolved_rulename_modifier
        self.local_rulename_modifier = local_rulename_modifier
        self.skip_validation = skip_validation
        self.skip_global_report_caption = skip_global_report_caption
        self.rule_whitelist = rule_whitelist
        self.rule_exclude_list = rule_exclude_list
        self.ruleinfo_overwrite = ruleinfo_overwrite
        self.allow_rule_overwrite = allow_rule_overwrite
        self.path_modifier = PathModifier(replace_prefix, prefix, workflow)
        self.replace_wrapper_tag = replace_wrapper_tag
        self.namespace = namespace
        self._ruleorder = _rules.Ruleorder()

    @property
    def rule_proxies(self) -> "_rules.Rules":
        return self.globals["rules"]

    def inherit_rule_proxies(self, child_modifier: "WorkflowModifier", stack_len: int):
        if self.rule_whitelist == []:
            for name, rule in child_modifier.rule_proxies._rules.items():
                if not rule._rescue:
                    self.rule_proxies._cache_rules[name] = (  # type: ignore[assignment]
                        self,
                        rule.rule.ruleinfo,
                        rule.rule.lineno,
                        rule.rule.snakefile,
                        rule.rule.is_checkpoint,
                        stack_len,
                    )
        else:
            for name, rule in child_modifier.rule_proxies._rules.items():
                name = child_modifier.local_rulename_modifier(name)
                self.rule_proxies._register_rule(name, rule)

    def inherit_ruleorder(
        self,
        child_modifier: "WorkflowModifier",
        resolved_rulename=_default_modify_rulename,
    ):
        for clause in child_modifier._ruleorder:
            _clause: list[str] = []
            for rulename in clause:
                modified_rulename = resolved_rulename(rulename)
                if modified_rulename in child_modifier.rule_proxies._rules:
                    _clause.append(modified_rulename)
            if len(_clause) > 1:
                self._ruleorder.add(*_clause)

    def skip_rule(self, rulename):
        return (
            self.rule_whitelist is not None and rulename not in self.rule_whitelist
        ) or (self.rule_exclude_list is not None and rulename in self.rule_exclude_list)

    def modify_rulename(self, rulename) -> str:
        return self.resolved_rulename_modifier(rulename)

    def modify_path(self, path, property=None):
        return self.path_modifier.modify(path, property)

    def modify_wrapper_uri(self, wrapper_uri, pattern=re.compile("^master/")):
        if self.replace_wrapper_tag is None or wrapper.is_url(wrapper_uri):
            return wrapper_uri
        else:
            return pattern.sub(self.replace_wrapper_tag + "/", wrapper_uri)

    def __enter__(self):
        # put this modifier on the stack, it becomes the currently valid modifier
        self.workflow.modifier_stack.append(self)

    def __exit__(self, type, value, traceback):
        # remove this modifier from the stack
        self.workflow.modifier_stack.pop()
