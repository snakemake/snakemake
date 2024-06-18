module deep_module:
    snakefile:
        "module_deep.smk"


use rule work from deep_module as deep_work


rule all:
    input:
        "foo.txt",
    default_target: True


# rules.deep_work has to work even if the rule is renamed in a parent module
# The rulename itself can be already modified.
assert hasattr(
    rules, "deep_work"
), f"bug: rule cannot be accessed as deep_work: {dir(rules)}"
