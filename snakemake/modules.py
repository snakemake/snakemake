import types

from snakemake.exceptions import WorkflowError


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
            skip_configfile=self.config is not None,
            skip_validation=self.skip_validation,
            rule_whitelist=self.get_rule_whitelist(rules),
            rulename_modifier=self.get_name_modifier_func(rules, name_modifier),
            ruleinfo_overwrite=ruleinfo,
            allow_rule_overwrite=True,
            namespace=self.name,
            replace_prefix=self.replace_prefix,
        ):
            self.workflow.include(snakefile, overwrite_first_rule=True)

    def get_name_modifier_func(self, rules=None, name_modifier=None):
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

    def get_snakefile(self):
        if self.meta_wrapper:
            # TODO add code to retrieve wrapper snakefile
            pass
        elif self.snakefile:
            return self.snakefile
        else:
            raise WorkflowError(
                "Module statement must either define snakefile or meta_wrapper to use."
            )

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
        skip_configfile=False,
        skip_validation=False,
        rulename_modifier=None,
        rule_whitelist=None,
        ruleinfo_overwrite=None,
        allow_rule_overwrite=False,
        replace_prefix=None,
        namespace=None,
    ):
        self.workflow = workflow

        self.globals = (
            globals if globals is not None else dict(workflow.vanilla_globals)
        )
        if config:
            self.globals["config"] = config

        self.skip_configfile = skip_configfile
        self.rulename_modifier = rulename_modifier
        self.skip_validation = skip_validation
        self.rule_whitelist = rule_whitelist
        self.ruleinfo_overwrite = ruleinfo_overwrite
        self.allow_rule_overwrite = allow_rule_overwrite
        self.path_modifier = PathModifier(replace_prefix) if replace_prefix else None
        self.namespace = namespace

    def skip_rule(self, rulename):
        return self.rule_whitelist is not None and rulename not in self.rule_whitelist

    def modify_rulename(self, rulename):
        if self.rulename_modifier is not None:
            return self.rulename_modifier(rulename)
        return rulename

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


class PathModifier:
    def __init__(self, replace_prefix: dict):
        self.skip_properties = set()

        import datrie
        self.trie = datrie.Trie("".join(set(char for prefix in replace_prefix for char in prefix)))
        for prefix, replacement in replace_prefix.items():
            self.trie[prefix] = replacement
    
    def modify(self, path: str, property: str):
        if property in self.skip_properties:
            return path
        prefixes = self.trie.prefix_items(path)
        if len(prefixes) > 1:
            # ambiguous prefixes
            raise WorkflowError("Multiple prefixes ({}) match the path {}. Make sure that the replace_prefix statement ""in your module definition does not yield ambiguous matches.".format(", ".join(prefix[0] for prefix in prefixes), path))
        elif prefixes:
            # replace prefix
            prefix, replacement = prefixes[0]
            return replacement + path[len(prefix):]
        else:
            # no matching prefix
            return path