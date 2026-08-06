__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
import types
import re
from typing import List, Optional, Set, Dict, Callable, TYPE_CHECKING
from snakemake.common import Rules

from snakemake.exceptions import CreateRuleException, WorkflowError
from snakemake.io.flags import DefaultFlags
from snakemake.path_modifier import PathModifier
from snakemake import wrapper
from snakemake.pathvars import Pathvars
from snakemake.logging import logger

if TYPE_CHECKING:
    from snakemake.workflow import Workflow
    from snakemake.sourcecache import SourceFile


def get_rule_whitelist(rules: list[str]):
    if "*" in rules:
        if len(rules) != 1:
            raise SyntaxError(
                "The 'use rule' statement uses a wildcard '*' but lists multiple additional rules."
            )
        return None
    return set(rules)


def get_name_modifier_func(rules: Optional[Set[str]], name_modifier=None):
    """
    Using for: use rule <rules> [from <module>] as <name_modifier>

    | rules     | name_modifier | effect
    | --------- | ------------- | ------
    | ["*"]     | None          | rule1, rule2, ...
    | ["*"]     | "A_*"         | A_rule1, A_rule2, ...
    | ["rule1"] | None          | rule1
    | ["rule1"] | "A_rule1"     | A_rule1
    | ["rule1"] | "A_*"         | A_rule1
    """
    if name_modifier is None:
        return None
    else:
        if "*" in name_modifier:
            return lambda rulename: name_modifier.replace("*", rulename)
        else:
            # Disallow constant renaming for wildcard or multi-rule imports.
            if not rules or len(rules) > 1:
                raise SyntaxError(
                    "Multiple rules in 'use rule' statement but name modification ('as' statement) does not contain a wildcard '*'."
                )
            return lambda rulename: name_modifier


class ModuleInfo:
    def __init__(
        self,
        workflow: "Workflow",
        name: str,
        pathvars: Pathvars,
        snakefile=None,
        meta_wrapper=None,
        config: "Optional[Dict]" = None,
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
        self.wildcards_modifier_overwrited: Dict[Optional[Callable], Set[str]] = {}
        self.pathvars = pathvars
        self.loaded = False

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

    @property
    def _cached_namespace(self) -> types.ModuleType:
        if not self.loaded:
            self.use_rules(set(), None)
            self.loaded = True
        return self.workflow.globals[self.name]

    @_cached_namespace.setter
    def _cached_namespace(self, value: types.ModuleType):
        if self.loaded:
            logger.warning(
                f"Reloading module {self.name}. "
                "This may cause problems if the module contains global variables"
                " or rules that are modified by the module. "
                "Consider using a unique name for each module import."
            )
        else:
            self.loaded = True
        self.workflow.globals[self.name] = value

    def use_rules(
        self,
        rules: Optional[Set],
        name_modifier: Optional[str] = None,
        exclude_rules: Optional[List[str]] = None,
        ruleinfo=None,
    ):
        """
        rules:
        - None -> load all rules,
        - empty set -> dry run (only load ruleinfo),
        - set of rule names -> (impossible)
            `Workflow.userule` performs rule-specific loading via `_cached_namespace`.
            so `ModuleInfo.use_rules` only sees `None` or empty set
        """
        self._cached_namespace = types.ModuleType(self.name)
        self._cached_namespace.__dict__.update(self.workflow.vanilla_globals)
        modifier = WorkflowModifier.for_module(
            self, rules, name_modifier, exclude_rules, ruleinfo
        )
        with modifier:
            self.workflow.include(
                modifier.base_snakefile, overwrite_default_target=True
            )
            if rules != set():
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

    def check_overwrite(self, rule_whitelist: Optional[Set[str]], rulename_modifier):
        """Check if the rule from module can be overwritten.

        > Specific rules may even be modified before using them,
        >  via a final with: followed by a block that lists items to overwrite.
        >  This modification can be performed after a general import,
        >  and will overwrite any unmodified import of the same rule.
        """
        if not rule_whitelist:  # all will be used, exclude rules do not matter
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
        workflow: "Workflow",
        globals: Dict,
        rule_proxies: Rules,
        path_modifier: PathModifier,
        pathvars: Pathvars,
        *,
        base_snakefile: "Optional[SourceFile]" = None,
        parent_modifier: "Optional[WorkflowModifier]" = None,
        resolved_rulename_modifier=None,
        ruleinfo_overwrite=None,
        allow_rule_overwrite=False,
    ):
        self.workflow = workflow
        self.globals = globals
        self.pathvars = pathvars
        self.rule_proxies = rule_proxies
        self.path_modifier = path_modifier

        self.base_snakefile = base_snakefile
        self.parent_modifier = parent_modifier
        self.resolved_rulename_modifier = resolved_rulename_modifier
        self.ruleinfo_overwrite = ruleinfo_overwrite
        self.allow_rule_overwrite = allow_rule_overwrite

        self.is_module = True
        self.namespace = None
        self.wildcard_constraints = dict()
        self.rules = set()
        self.modules: Dict[str, ModuleInfo] = dict()

        self.skip_configfile = False
        self.skip_validation = False
        self.skip_global_report_caption = False
        self.rule_whitelist: "Optional[Set[str]]" = None
        self.rule_exclude_list: "Optional[List[str]]" = None
        self.replace_wrapper_tag = None

        self.default_input_flags: DefaultFlags = DefaultFlags()
        self.default_output_flags: DefaultFlags = DefaultFlags()

    def _post_set_globals(self, config=None):
        self.globals["__name__"] = self.namespace
        self.globals["config"] = {} if config is None else config
        self.globals["rules"] = self.rule_proxies
        self.globals["checkpoints"] = self.globals["checkpoints"].spawn_new_namespace()

    @classmethod
    def for_modifier0(cls, workflow, globals: Dict):
        self = cls(
            workflow,
            globals=globals,
            rule_proxies=Rules(),
            path_modifier=PathModifier(None, None, workflow),
            pathvars=Pathvars.with_defaults(),
        )
        self._post_set_globals()
        return self

    @classmethod
    def for_module(
        cls,
        module_info: ModuleInfo,
        rules: Optional[Set],
        name_modifier: "Optional[str]" = None,
        exclude_rules: "Optional[List[str]]" = None,
        ruleinfo=None,
    ):
        rulename_modifier = get_name_modifier_func(rules, name_modifier)
        workflow = module_info.workflow
        module_info.wildcards_modifier_overwrited[rulename_modifier] = set()
        self = cls(
            workflow,
            globals=module_info._cached_namespace.__dict__,
            rule_proxies=module_info.rule_proxies,
            path_modifier=module_info.path_modifier,
            pathvars=module_info.pathvars,
            base_snakefile=module_info.get_snakefile(),
            parent_modifier=workflow.modifier,
            resolved_rulename_modifier=rulename_modifier,
            ruleinfo_overwrite=ruleinfo,
            allow_rule_overwrite=workflow.modifier.allow_rule_overwrite,
        )
        self.namespace = module_info.name
        self._post_set_globals(module_info.config)
        self.skip_configfile = module_info.config is not None
        self.skip_validation = module_info.skip_validation
        self.skip_global_report_caption = (
            workflow.report_text is not None
        )  # do not overwrite existing report text via module
        self.rule_whitelist = rules
        self.rule_exclude_list = exclude_rules
        self.replace_wrapper_tag = module_info.get_wrapper_tag()
        return self

    @classmethod
    def for_userule(
        cls,
        workflow,
        globals: Dict,
        pathvars: Pathvars,
        snakefile,
        ruleinfo=None,
        allow_overwrite=False,
    ):
        parent_modifier = workflow.modifier
        self = cls(
            workflow,
            globals=globals,
            rule_proxies=parent_modifier.rule_proxies,
            path_modifier=parent_modifier.path_modifier,
            pathvars=pathvars,
            base_snakefile=snakefile,
            parent_modifier=parent_modifier,
            resolved_rulename_modifier=None,
            ruleinfo_overwrite=ruleinfo,
            allow_rule_overwrite=allow_overwrite | parent_modifier.allow_rule_overwrite,
        )
        self.is_module = False
        self.wildcard_constraints = parent_modifier.wildcard_constraints
        self.rules = parent_modifier.rules
        self.modules = parent_modifier.modules
        return self

    def is_main_snakefile(self):
        return self.base_snakefile is None

    def _skip_rulename(self, rulename):
        if self.rule_whitelist is not None and rulename not in self.rule_whitelist:
            return True
        if self.rule_exclude_list is not None and rulename in self.rule_exclude_list:
            return True
        return False

    def inherit_rule_proxies(self, child_modifier: "WorkflowModifier"):
        for name, rule in child_modifier.rule_proxies._rules.items():
            if child_modifier._skip_rulename(name):
                continue
            name1 = child_modifier._modify_rulename(name)
            self.rule_proxies._register_rule(name1, rule)

    def avail_rulename(self, rulename):
        if self._skip_rulename(rulename):
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
        if not self.is_module and self.base_snakefile is not None:
            self.workflow.included_stack.append(self.base_snakefile)

    def __exit__(self, exc_type, exc_value, traceback):
        # remove this modifier from the stack
        self.workflow.modifier_stack.pop()
        if not self.is_module and self.base_snakefile is not None:
            self.workflow.included_stack.pop()
