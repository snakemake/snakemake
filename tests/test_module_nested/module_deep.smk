a = 1


rule run:
    output:
        txt=touch("A/foo.txt"),
    shell:
        f"echo 'I was here' > {output.txt}"


rule run2:
    input:
        rds=rules.run.output.txt,
    output:
        txt="output/run2.txt",
    shell:
        "cp {input.rds} {output.txt}"


use rule run as run3 with:
    output:
        "0",


# rules.run has to run even if the rule is renamed in a parent module
# The rulename itself can be already modified.
assert hasattr(rules, "run"), f"bug: rule cannot be accessed as run: {dir(rules)}"
