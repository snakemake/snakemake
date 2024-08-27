__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from collections import namedtuple
from copy import copy
from pathlib import Path

from .settings.types import ResourceSettings, WorkflowSettings
from .logging import logger
from .rules import Rule
from .exceptions import RuleException
from .resources import ParsedResource
from .wrapper import get_conda_env
from . import modules

InOutput = namedtuple("InOutput", ["paths", "kwpaths", "modifier"])


def get_resource_value(value):
    if isinstance(value, ParsedResource):
        return value.value
    else:
        return value


class RuleInfo:
    ref_attributes = {"func", "path_modifier"}

    def __init__(self, func=None):
        self.func = func
        self.shellcmd = None
        self.name = None
        self.norun = False
        self.input = None
        self.output = None
        self.params = None
        self.message = None
        self.benchmark = None
        self.conda_env = None
        self.container_img = None
        self.is_containerized = False
        self.env_modules = None
        self.wildcard_constraints = None
        self.threads = None
        self.shadow_depth = None
        self.resources = None
        self.priority = None
        self.retries = None
        self.log = None
        self.docstring = None
        self.group = None
        self.script = None
        self.notebook = None
        self.wrapper = None
        self.template_engine = None
        self.cwl = None
        self.cache = False
        self.path_modifier = None
        self.handover = False
        self.default_target = False
        self.localrule = False

    def __copy__(self):
        """Return a copy of this ruleinfo."""
        ruleinfo = RuleInfo(self.func)
        for attribute in self.__dict__:
            if attribute in self.ref_attributes:
                setattr(ruleinfo, attribute, getattr(self, attribute))
            else:
                # shallow copies are enough
                setattr(ruleinfo, attribute, copy(getattr(self, attribute)))
        return ruleinfo

    def apply_modifier(
        self,
        modifier: "modules.WorkflowModifier",
        prefix_replacables={"input", "output", "log", "benchmark"},
    ):
        """Update this ruleinfo with the given one (used for 'use rule' overrides)."""
        path_modifier = modifier.path_modifier
        skips = set()

        if modifier.ruleinfo_overwrite:
            for key, value in modifier.ruleinfo_overwrite.__dict__.items():
                if key != "func" and value is not None:
                    self.__dict__[key] = value
                    if key in prefix_replacables:
                        skips.add(key)

        if path_modifier.modifies_prefixes and skips:
            # use a specialized copy of the path modifier
            path_modifier = copy(path_modifier)
            path_modifier.skip_properties = skips
        # add path modifier
        self.path_modifier = path_modifier

        # modify wrapper if requested
        self.wrapper = modifier.modify_wrapper_uri(self.wrapper)

    @property
    def is_invalid(ruleinfo):
        return not (
            ruleinfo.script
            or ruleinfo.wrapper
            or ruleinfo.shellcmd
            or ruleinfo.notebook
        )

    def update_rule(ruleinfo, rule: Rule):
        if ruleinfo.wildcard_constraints:
            rule.set_wildcard_constraints(
                *ruleinfo.wildcard_constraints[0],
                **ruleinfo.wildcard_constraints[1],
            )
        if ruleinfo.input:
            rule.input_modifier = ruleinfo.input.modifier
            rule.set_input(*ruleinfo.input.paths, **ruleinfo.input.kwpaths)
        if ruleinfo.output:
            rule.output_modifier = ruleinfo.output.modifier
            rule.set_output(*ruleinfo.output.paths, **ruleinfo.output.kwpaths)
        if ruleinfo.log:
            rule.log_modifier = ruleinfo.log.modifier
            rule.set_log(*ruleinfo.log.paths, **ruleinfo.log.kwpaths)
        if ruleinfo.params:
            rule.set_params(*ruleinfo.params[0], **ruleinfo.params[1])

        if ruleinfo.shadow_depth:
            if ruleinfo.shadow_depth not in (
                True,
                "shallow",
                "full",
                "minimal",
                "copy-minimal",
            ):
                raise RuleException(
                    "Shadow must either be 'minimal', 'copy-minimal', 'shallow', 'full', "
                    "or True (equivalent to 'full')",
                    rule=rule,
                )
            if ruleinfo.shadow_depth is True:
                rule.shadow_depth = "full"
                logger.warning(
                    f"Shadow is set to True in rule {rule} (equivalent to 'full'). "
                    "It's encouraged to use the more explicit options "
                    "'minimal|copy-minimal|shallow|full' instead."
                )
            else:
                rule.shadow_depth = ruleinfo.shadow_depth

        if ruleinfo.resources:
            resources: dict
            args, resources = ruleinfo.resources
            if args:
                raise RuleException("Resources have to be named.")
            if not all(
                map(
                    lambda r: isinstance(r, int) or isinstance(r, str) or callable(r),
                    resources.values(),
                )
            ):
                raise RuleException(
                    "Resources values have to be integers, strings, or callables (functions)",
                    rule=rule,
                )
            rule.resources.update(resources)

        if ruleinfo.priority:
            if not isinstance(ruleinfo.priority, int) and not isinstance(
                ruleinfo.priority, float
            ):
                raise RuleException("Priority values have to be numeric.", rule=rule)
            rule.priority = ruleinfo.priority

        if ruleinfo.retries:
            if not isinstance(ruleinfo.retries, int) or ruleinfo.retries < 0:
                raise RuleException(
                    "Retries values have to be integers >= 0", rule=rule
                )
        rule.restart_times = ruleinfo.retries

        if ruleinfo.message:
            rule.message = ruleinfo.message
        if ruleinfo.benchmark:
            rule.benchmark_modifier = ruleinfo.benchmark.modifier
            rule.benchmark = ruleinfo.benchmark.paths

        group = ruleinfo.group
        if group is not None:
            rule.group = group

        if ruleinfo.env_modules:
            # If using environment modules and they are defined for the rule,
            # ignore conda and singularity directive below.
            # The reason is that this is likely intended in order to use
            # a software stack specifically compiled for a particular
            # HPC cluster.
            if ruleinfo.is_invalid:
                raise RuleException(
                    "envmodules directive is only allowed with "
                    "shell, script, notebook, or wrapper directives (not with run or the template_engine)",
                    rule=rule,
                )
            from .deployment.env_modules import EnvModules

            rule.env_modules = EnvModules(*ruleinfo.env_modules)
        if ruleinfo.conda_env:
            if ruleinfo.is_invalid:
                raise RuleException(
                    "Conda environments are only allowed "
                    "with shell, script, notebook, or wrapper directives "
                    "(not with run or template_engine).",
                    rule=rule,
                )
            if isinstance(ruleinfo.conda_env, Path):
                ruleinfo.conda_env = str(ruleinfo.conda_env)
            rule.conda_env = ruleinfo.conda_env

        if ruleinfo.container_img:
            if ruleinfo.is_invalid:
                raise RuleException(
                    "Singularity directive is only allowed "
                    "with shell, script, notebook or wrapper directives "
                    "(not with run or template_engine).",
                    rule=rule,
                )
            rule.container_img = ruleinfo.container_img
            rule.is_containerized = ruleinfo.is_containerized

        rule.norun = ruleinfo.norun
        rule.docstring = ruleinfo.docstring
        rule.run_func = ruleinfo.func
        rule.shellcmd = ruleinfo.shellcmd
        rule.script = ruleinfo.script
        rule.notebook = ruleinfo.notebook
        rule.wrapper = ruleinfo.wrapper
        rule.template_engine = ruleinfo.template_engine
        rule.cwl = ruleinfo.cwl
        rule.ruleinfo = ruleinfo

    def update_rule_settings(
        ruleinfo,
        rule: Rule,
        resource_settings=ResourceSettings(),
        workflow_settings=WorkflowSettings(),
        container_img=None,
        is_containerized=False,
    ):
        name = rule.name
        if ruleinfo.threads is not None:
            if (
                not isinstance(ruleinfo.threads, int)
                and not isinstance(ruleinfo.threads, float)
                and not callable(ruleinfo.threads)
            ):
                raise RuleException(
                    "Threads value has to be an integer, float, or a callable.",
                    rule=rule,
                )
            if name not in resource_settings.overwrite_threads:
                if isinstance(ruleinfo.threads, float):
                    ruleinfo.threads = int(ruleinfo.threads)
                rule.resources["_cores"] = ruleinfo.threads
        else:
            rule.resources["_cores"] = 1

        if name in resource_settings.overwrite_threads:
            rule.resources["_cores"] = get_resource_value(
                resource_settings.overwrite_threads[name]
            )
        if name in resource_settings.overwrite_resources:
            rule.resources.update(
                (resource, get_resource_value(value))
                for resource, value in resource_settings.overwrite_resources[
                    name
                ].items()
            )

        if ruleinfo.wrapper:
            rule.conda_env = get_conda_env(
                ruleinfo.wrapper, prefix=workflow_settings.wrapper_prefix
            )
            # TODO retrieve suitable singularity image

        if not ruleinfo.container_img and container_img:
            if not ruleinfo.is_invalid and ruleinfo.container_img != False:
                # skip rules with run directive or empty image
                rule.container_img = container_img
                rule.is_containerized = is_containerized
