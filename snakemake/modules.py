__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from contextlib import contextmanager
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
        with self.modifier.mask_rules():
            self.workflow.include(snakefile, overwrite_default_target=False)
            self.stack_len = len(self.workflow.included_stack) + 1
        if self.name:
            # update globals in inner modules
            self.namespace.__dict__.update(self.modifier.globals)
            self.workflow.globals[self.modifier.namespace] = self.namespace
            for _modifier, *_ in self.rule_proxies._cached.values():
                # all rules in _cache now are defined in the module,
                #  or used from the submodules.
                # They cannot be skipped.
                _modifier.rule_whitelist = None

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
        avail_rules = self.get_rules(rules, exclude_rules)
        with self.modifier.mask_ruleinfo(ruleinfo):
            self._include(
                avail_rules,
                self.snakefile,
                self.rule_proxies._stack,
                name_modifier_func,
                True,
            )
        self.parent_modifier.inherit_rule_proxies(self.modifier, self.stack_len)
        self.parent_modifier.inherit_ruleorder(self.modifier, name_modifier_func)

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

    def get_rules(
        self,
        rule_whitelist: list[str],
        exclude_rules: list[str] | None = None,
    ):
        """
        Here, rules are reordered back as loaded from the module.
        It is used to inherent correct default rule from the module.
        """
        excludes = set() if exclude_rules is None else set(exclude_rules)
        if "*" in rule_whitelist:
            if len(rule_whitelist) != 1:
                raise SyntaxError(
                    "The 'use rule' statement uses a wildcard '*' but lists multiple additional rules."
                )
            return self.rule_proxies._cached.keys() - excludes
        return set(rule_whitelist) - excludes

    def _include(
        self,
        avail_rules: set[str],
        snakefile: str | None,
        rule_stacks: list[str | tuple[str, list]],
        name_modifier_func=_default_modify_rulename,
        overwrite_default_target=False,
    ):
        """
        This include is designed for select the default_target
         as the old version does.
        The param `overwrite_default_target` should only be assigned as True for the base level

        Each rule will be send to Workflow with it's own modifier
        """
        with self.workflow._include_stack(
            snakefile,
            overwrite_default_target=overwrite_default_target,
            module_use=True,
        ):
            for rule_s in rule_stacks:
                if isinstance(rule_s, tuple):
                    self._include(avail_rules, rule_s[0], rule_s[1], name_modifier_func)
                elif rule_s in avail_rules:
                    resolved_rulename = name_modifier_func(rule_s)
                    modifier, ruleinfo_, lineno, snakefile, checkpoint, stack_len = (
                        self.rule_proxies._cached[rule_s]
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
        rule_whitelist: list | None = None,
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
            """
            During module loading, may use rule from submodules.
            Then we just need to add the "used" rule to cache only
            """
            for name, rule in child_modifier.rule_proxies._used.items():
                self.rule_proxies._cached[name] = (  # type: ignore[assignment]
                    self,
                    rule.rule.ruleinfo,
                    rule.rule.lineno,
                    rule.rule.snakefile,
                    rule.rule.is_checkpoint,
                    stack_len,
                )
            self.rule_proxies._stack.append(
                (child_modifier.base_snakefile, list(child_modifier.rule_proxies._used))
            )
        else:
            for name, rule in child_modifier.rule_proxies._used.items():
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
                if modified_rulename in child_modifier.rule_proxies._used:
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

    @contextmanager
    def mask_rules(self):
        workflow_rules, self.workflow._rules = self.workflow._rules, {}
        try:
            with self:
                yield
        finally:
            self.workflow._rules = workflow_rules

    @contextmanager
    def mask_skip(self):
        white, self.rule_whitelist = self.rule_whitelist, None
        exclude, self.rule_exclude_list = (
            self.rule_exclude_list,
            None,
        )
        try:
            with self:
                yield
        finally:
            self.rule_whitelist = white
            self.rule_exclude_list = exclude

    @contextmanager
    def mask_ruleinfo(self, ruleinfo: "_workflow.RuleInfo | None" = None):
        ruleinfo, self.ruleinfo_overwrite = self.ruleinfo_overwrite, None
        try:
            with self:
                yield
        finally:
            self.ruleinfo_overwrite = ruleinfo
