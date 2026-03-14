shell.executable("bash")


module A:
    snakefile:
        "module_deep.smk"
    prefix:
        "b-"


use rule * from A as A_*


if "bb" not in config:

    use rule run from A as A_run with:
        threads: 1


use rule run from A as A_run with:
    input:
        rds=config.get("bb", ".snakemake"),
    log:
        "aaalog",


rule all:
    input:
        "foo.txt",
    default_target: True


# rules.deep_work has to work even if the rule is renamed in a parent module
# The rulename itself can be already modified.
assert hasattr(
    rules, "A_run"
), f"bug: rule cannot be accessed as A_run: {dir(rules)}"
