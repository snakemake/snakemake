__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from collections import namedtuple
from copy import copy


InOutput = namedtuple("InOutput", ["paths", "kwpaths", "modifier"])


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
        self, modifier, prefix_replacables={"input", "output", "log", "benchmark"}
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
