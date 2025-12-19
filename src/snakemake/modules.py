__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
import types
import re
from typing import List, Optional, Set, Dict, Callable
from snakemake.common import Rules

from snakemake.exceptions import CreateRuleException, WorkflowError
from snakemake.io.flags import DefaultFlags
from snakemake.path_modifier import PathModifier
from snakemake import wrapper
from snakemake.pathvars import Pathvars


def get_name_modifier_func(rules: List[str], name_modifier=None):
    if name_modifier is None:
        return None
    else:
        if "*" in name_modifier:
            return lambda rulename: name_modifier.replace("*", rulename)
        elif name_modifier is not None:
            # Disallow constant renaming for wildcard or multi-rule imports.
            if "*" in rules or len(rules) > 1:
                raise SyntaxError(
                    "Multiple rules in 'use rule' statement but name modification ('as' statement) does not contain a wildcard '*'."
                )
            return lambda rulename: name_modifier


class ModuleInfo:
    def __init__(
        self,
        workflow,
        name,
        pathvars: Pathvars,
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
        self.parent_modifier: WorkflowModifier = self.workflow.modifier
        self.rule_proxies = Rules()
        self.wildcards_modifier_overwrited: Dict[Callable | None, Set[str]] = {}
        self.pathvars = pathvars

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
        self.path_modifier = PathModifier(
            replace_prefix, prefix, workflow, self.parent_modifier.path_modifier
        )

    def use_rules(
        self,
        rules: List[str],
        name_modifier: str | None = None,
        exclude_rules: List[str] | None = None,
        ruleinfo=None,
        skip_global_report_caption=False,
    ):
        snakefile = self.get_snakefile()
        rule_whitelist = self.get_rule_whitelist(rules)
        rulename_modifier = get_name_modifier_func(rules, name_modifier)
        modifier = WorkflowModifier(
            self.workflow,
            config=self.config,
            base_snakefile=snakefile,
            skip_validation=self.skip_validation,
            skip_global_report_caption=skip_global_report_caption,
            rule_exclude_list=exclude_rules,
            rule_whitelist=rule_whitelist,
            resolved_rulename_modifier=rulename_modifier,
            local_rulename_modifier=rulename_modifier,
            ruleinfo_overwrite=ruleinfo,
            allow_rule_overwrite=self.check_overwrite(
                rule_whitelist, rulename_modifier
            ),
            namespace=self.name,
            path_modifier=self.path_modifier,
            replace_wrapper_tag=self.get_wrapper_tag(),
            rule_proxies=self.rule_proxies,
            pathvars=self.pathvars,
        )
        with modifier:
            self.workflow.include(snakefile, overwrite_default_target=True)
            self.parent_modifier.inherit_rule_proxies(modifier)

    def get_snakefile(self):
        if self.meta_wrapper:
            for snakefile in ("/meta_wrapper.smk", "/test/Snakefile"):
                path = wrapper.get_path(
                    self.meta_wrapper + snakefile,
                    self.workflow.workflow_settings.wrapper_prefix,
                )
                if self.workflow.sourcecache.exists(path):
                    return path
            raise WorkflowError(
                f"Invalid meta wrapper {self.meta_wrapper}: Could not find "
                "meta_wrapper.smk or test/Snakefile (old style)."
            )
        elif self.snakefile:
            return self.snakefile
        else:
            raise WorkflowError(
                "Module statement must either define snakefile or meta_wrapper to use."
            )

    def get_wrapper_tag(self):
        from packaging.version import Version

        if self.meta_wrapper:
            if wrapper.is_url(self.meta_wrapper):
                # no wrapper tag replacement, use meta-wrapper as is
                return None
            tag = self.meta_wrapper.split("/", 1)[0]
            ver_match = wrapper.ver_regex.match(tag)
            if ver_match and Version(ver_match.group("ver")) >= Version("8.0.0"):
                # New style meta-wrappers, containing concrete versions of each wrapper
                return None
            else:
                return tag
        return None

    def get_rule_whitelist(self, rules):
        if "*" in rules:
            if len(rules) != 1:
                raise SyntaxError(
                    "The 'use rule' statement uses a wildcard '*' but lists multiple additional rules."
                )
            return None
        return set(rules)

    def check_overwrite(self, rule_whitelist: Set[str] | None, rulename_modifier):
        """Check if the rule from module can be overwritten.

        > Specific rules may even be modified before using them,
        >  via a final with: followed by a block that lists items to overwrite.
        >  This modification can be performed after a general import,
        >  and will overwrite any unmodified import of the same rule.
        """
        if not rule_whitelist:  # all will be used, exclude rules do not matter
            self.wildcards_modifier_overwrited[rulename_modifier] = set()
            return False
        matched = set()
        for rule in rule_whitelist:
            modified_this = rulename_modifier(rule) if rulename_modifier else rule
            for modified in self.wildcards_modifier_overwrited:
                if (modified(rule) if modified else rule) == modified_this:
                    if rule in self.wildcards_modifier_overwrited[modified]:
                        raise CreateRuleException(
                            f"The rule {rule} is imported with same name modifier and modified more than once",
                            snakefile=self.snakefile,
                        )
                    matched.add(rule)
                    self.wildcards_modifier_overwrited[modified].add(rule)
        # Now, there is already a rule with the same name defined before.
        # We confirm that it is created by `use rule * ...`
        # So we can safely overwrite it.
        return matched == set(rule_whitelist)


class WorkflowModifier:
    def __init__(
        self,
        workflow,
        is_module=True,
        globals=None,
        config=None,
        base_snakefile=None,
        skip_validation=False,
        skip_global_report_caption=False,
        resolved_rulename_modifier=None,
        local_rulename_modifier=None,
        rule_whitelist: Set[str] | None = None,
        rule_exclude_list=None,
        ruleinfo_overwrite=None,
        allow_rule_overwrite=False,
        path_modifier: PathModifier | None = None,
        replace_wrapper_tag=None,
        namespace=None,
        rule_proxies: Rules | None = None,
        pathvars: Optional[Pathvars] = None,
    ):
        if is_module:
            if globals is None:  # use rule from module with maybe_ruleinfo
                globals = dict(workflow.vanilla_globals)
                self.parent_modifier = workflow.modifier
                allow_rule_overwrite |= self.parent_modifier.allow_rule_overwrite
            else:
                # the first module modifier of workflow
                self.parent_modifier = None
            self.globals = globals
            self.globals["__name__"] = namespace
            self.globals["rules"] = self.rule_proxies = rule_proxies or Rules()
            self.globals["checkpoints"] = self.globals[
                "checkpoints"
            ].spawn_new_namespace()

            self.globals["config"] = config if config is not None else {}
            self.wildcard_constraints: dict = dict()

            assert (
                pathvars is not None
            )  # pathvars is only None in case of is_module=False
            self.pathvars = pathvars
            self.rules: set = set()
            self.modules: dict = dict()
            self.path_modifier = path_modifier or PathModifier(None, None, workflow)
        else:
            # use rule (from same include) as ... with: init with values from parent modifier
            self.parent_modifier = parent_modifier = workflow.modifier
            self.globals = parent_modifier.globals
            self.wildcard_constraints = parent_modifier.wildcard_constraints
            self.pathvars = parent_modifier.pathvars
            self.rules = parent_modifier.rules
            self.rule_proxies = parent_modifier.rule_proxies
            self.modules = parent_modifier.modules
            self.path_modifier = parent_modifier.path_modifier
            allow_rule_overwrite |= self.parent_modifier.allow_rule_overwrite

        self.is_module = is_module
        self.workflow = workflow
        self.base_snakefile = base_snakefile

        self.skip_configfile = config is not None
        self.resolved_rulename_modifier = resolved_rulename_modifier
        self.local_rulename_modifier = local_rulename_modifier
        self.skip_validation = skip_validation
        self.skip_global_report_caption = skip_global_report_caption
        self.rule_whitelist = rule_whitelist
        self.rule_exclude_list = rule_exclude_list
        self.ruleinfo_overwrite = ruleinfo_overwrite
        self.allow_rule_overwrite = allow_rule_overwrite
        self.replace_wrapper_tag = replace_wrapper_tag
        self.namespace = namespace
        self.default_input_flags: DefaultFlags = DefaultFlags()
        self.default_output_flags: DefaultFlags = DefaultFlags()

    def inherit_rule_proxies(self, child_modifier: "WorkflowModifier"):
        for name, rule in child_modifier.rule_proxies._rules.items():
            if child_modifier.local_rulename_modifier is not None:
                name = child_modifier.local_rulename_modifier(name)
            self.rule_proxies._register_rule(name, rule)

    def avail_rulename(self, rulename):
        if (
            self.rule_whitelist is not None and rulename not in self.rule_whitelist
        ) or (
            self.rule_exclude_list is not None and rulename in self.rule_exclude_list
        ):
            return None
        if self.parent_modifier is not None:
            name_modifier = self._modify_rulename(rulename)
            return self.parent_modifier.avail_rulename(name_modifier)
        return rulename

    def modify_rulename(self, rulename):
        rulename = self._modify_rulename(rulename)
        if self.parent_modifier is not None:
            rulename = self.parent_modifier.modify_rulename(rulename)
        return rulename

    def _modify_rulename(self, rulename):
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
