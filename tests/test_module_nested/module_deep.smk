# deep_module.smk
rule all:
    input:
        "foo.txt",


rule work:
    output:
        "foo.txt",
    shell:
        "echo 'I was here' > {output}"


assert rules.work.rule.name == "work", rules.work.rule.name
