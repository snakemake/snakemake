# deep_module.smk
rule all:
    input:
        "foo.txt",


rule work:
    output:
        "foo.txt",
    shell:
        "echo 'I was here' > {output}"


# rules.work has to work even if the rule is renamed in a parent module
# The rulename itself can be already modified.
assert hasattr(rules, "work"), f"bug: rule cannot be accessed as work: {dir(rules)}"
