TABLES = "results/tables/"

rule a:
    input:
        "test.in"
    output:
        TABLES + "test.out"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"