__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
import types
import re
from snakemake.common import Rules

from snakemake.exceptions import WorkflowError
from snakemake.path_modifier import PathModifier
from snakemake import wrapper


def get_name_modifier_func(rules=None, name_modifier=None, parent_modifier=None):
    if name_modifier is None:
        return None
    else:
        if parent_modifier is None:
            parent_modifier_func = lambda rulename: rulename
        else:
            parent_modifier_func = parent_modifier.modify_rulename
        if "*" in name_modifier:
            return lambda rulename: parent_modifier_func(
                name_modifier.replace("*", rulename)
            )
        elif name_modifier is not None:
            if len(rules) > 1:
                raise SyntaxError(
                    "Multiple rules in 'use rule' statement but name modification ('as' statement) does not contain a wildcard '*'."
                )
            return lambda rulename: parent_modifier_func(name_modifier)


class ModuleInfo:
    def __init__(
        self,
        workflow,
        name,
        snakefile=None,
        meta_wrapper=None,
        config=None,
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
        self.rule_proxies = Rules()

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

    def use_rules(
        self,
        rules=None,
        name_modifier=None,
        exclude_rules=None,
        ruleinfo=None,
        skip_global_report_caption=False,
    ):
        snakefile = self.get_snakefile()
        modifier = WorkflowModifier(
            self.workflow,
            config=self.config,
            base_snakefile=snakefile,
            skip_configfile=self.config is not None,
            skip_validation=self.skip_validation,
            skip_global_report_caption=skip_global_report_caption,
            rule_exclude_list=exclude_rules,
            rule_whitelist=self.get_rule_whitelist(rules),
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
        with modifier:
            self.workflow.include(snakefile, overwrite_default_target=True)
            self.parent_modifier.inherit_rule_proxies(modifier)

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

    def get_rule_whitelist(self, rules):
        if "*" in rules:
            if len(rules) != 1:
                raise SyntaxError(
                    "The 'use rule' statement uses a wildcard '*' but lists multiple additional rules."
                )
            else:
                return None
        return set(rules)


class WorkflowModifier:
    def __init__(
        self,
        workflow,
        parent_modifier=None,
        globals=None,
        config=None,
        base_snakefile=None,
        skip_configfile=False,
        skip_validation=False,
        skip_global_report_caption=False,
        resolved_rulename_modifier=None,
        local_rulename_modifier=None,
        rule_whitelist=None,
        rule_exclude_list=None,
        ruleinfo_overwrite=None,
        allow_rule_overwrite=False,
        replace_prefix=None,
        prefix=None,
        replace_wrapper_tag=None,
        namespace=None,
        rule_proxies=None,
    ):
        if parent_modifier is not None:
            # init with values from parent modifier
            self.base_snakefile = parent_modifier.base_snakefile
            self.globals = parent_modifier.globals
            self.skip_configfile = parent_modifier.skip_configfile
            self.resolved_rulename_modifier = parent_modifier.resolved_rulename_modifier
            self.local_rulename_modifier = parent_modifier.local_rulename_modifier
            self.skip_validation = parent_modifier.skip_validation
            self.skip_global_report_caption = parent_modifier.skip_global_report_caption
            self.rule_whitelist = parent_modifier.rule_whitelist
            self.rule_exclude_list = parent_modifier.rule_exclude_list
            self.ruleinfo_overwrite = parent_modifier.ruleinfo_overwrite
            self.allow_rule_overwrite = parent_modifier.allow_rule_overwrite
            self.path_modifier = parent_modifier.path_modifier
            self.replace_wrapper_tag = parent_modifier.replace_wrapper_tag
            self.namespace = parent_modifier.namespace
            self.wildcard_constraints = parent_modifier.wildcard_constraints
            self.rules = parent_modifier.rules
            self.rule_proxies = parent_modifier.rule_proxies
        else:
            # default settings for globals if not inheriting from parent
            self.globals = (
                globals if globals is not None else dict(workflow.vanilla_globals)
            )
            self.wildcard_constraints = dict()
            self.rules = set()
            self.rule_proxies = rule_proxies or Rules()
            self.globals["rules"] = self.rule_proxies

        self.workflow = workflow
        self.base_snakefile = base_snakefile

        if config is not None:
            self.globals["config"] = config

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

    def inherit_rule_proxies(self, child_modifier):
        for name, rule in child_modifier.rule_proxies._rules.items():
            if child_modifier.local_rulename_modifier is not None:
                name = child_modifier.local_rulename_modifier(name)
            self.rule_proxies._register_rule(name, rule)

    def skip_rule(self, rulename):
        return (
            self.rule_whitelist is not None and rulename not in self.rule_whitelist
        ) or (self.rule_exclude_list is not None and rulename in self.rule_exclude_list)

    def modify_rulename(self, rulename):
        if self.resolved_rulename_modifier is not None:
            return self.resolved_rulename_modifier(rulename)
        return rulename

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
        if self.namespace:
            namespace = types.ModuleType(self.namespace)
            namespace.__dict__.update(self.globals)
            self.workflow.globals[self.namespace] = namespace
