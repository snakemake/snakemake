module deep_module:
    snakefile:
        "module_deep.smk"


use rule work from deep_module as deep_work


rule all:
    input:
        "foo.txt",
    default_target: True


assert rules.deep_work.rule.name == "deep_work", rules.deep_work.rule.name
