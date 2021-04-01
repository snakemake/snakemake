__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import types
import re

from snakemake.exceptions import WorkflowError
from snakemake.path_modifier import PathModifier
from snakemake import wrapper
from snakemake.checkpoints import Checkpoints
from snakemake.common import Rules, Scatter, Gather


def get_name_modifier_func(rules=None, name_modifier=None):
    if name_modifier is None:
        return None
    else:
        if "*" in name_modifier:
            return lambda rulename: name_modifier.replace("*", rulename)
        elif name_modifier is not None:
            if len(rules) > 1:
                raise SyntaxError(
                    "Multiple rules in 'use rule' statement but name modification ('as' statement) does not contain a wildcard '*'."
                )
            return lambda rulename: name_modifier


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
    ):
        self.workflow = workflow
        self.name = name
        self.snakefile = snakefile
        self.meta_wrapper = meta_wrapper
        self.config = config
        self.skip_validation = skip_validation
        self.replace_prefix = replace_prefix

    def use_rules(self, rules=None, name_modifier=None, ruleinfo=None):
        snakefile = self.get_snakefile()
        with WorkflowModifier(
            self.workflow,
            config=self.config,
            base_snakefile=snakefile,
            skip_configfile=self.config is not None,
            skip_validation=self.skip_validation,
            rule_whitelist=self.get_rule_whitelist(rules),
            rulename_modifier=get_name_modifier_func(rules, name_modifier),
            ruleinfo_overwrite=ruleinfo,
            allow_rule_overwrite=True,
            namespace=self.name,
            replace_prefix=self.replace_prefix,
            replace_wrapper_tag=self.get_wrapper_tag(),
        ):
            self.workflow.include(snakefile, overwrite_first_rule=True)

    def get_snakefile(self):
        if self.meta_wrapper:
            return wrapper.get_path(
                self.meta_wrapper + "/test/Snakefile", self.workflow.wrapper_prefix
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
        globals=None,
        config=None,
        base_snakefile=None,
        skip_configfile=False,
        skip_validation=False,
        rulename_modifier=None,
        rule_whitelist=None,
        ruleinfo_overwrite=None,
        allow_rule_overwrite=False,
        replace_prefix=None,
        replace_wrapper_tag=None,
        namespace=None,
    ):
        self.workflow = workflow
        self.base_snakefile = base_snakefile

        self.globals = (
            globals if globals is not None else dict(workflow.vanilla_globals)
        )
        if config is not None:
            self.globals["config"] = config

        self.skip_configfile = skip_configfile
        self.rulename_modifier = rulename_modifier
        self.skip_validation = skip_validation
        self.rule_whitelist = rule_whitelist
        self.ruleinfo_overwrite = ruleinfo_overwrite
        self.allow_rule_overwrite = allow_rule_overwrite
        self.path_modifier = PathModifier(replace_prefix, workflow)
        self.replace_wrapper_tag = replace_wrapper_tag
        self.namespace = namespace

    def skip_rule(self, rulename):
        return self.rule_whitelist is not None and rulename not in self.rule_whitelist

    def modify_rulename(self, rulename):
        if self.rulename_modifier is not None:
            return self.rulename_modifier(rulename)
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
