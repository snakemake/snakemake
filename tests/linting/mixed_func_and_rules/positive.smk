include: "common.smk"

rule a:
    input:
        get_dataset
    output:
        "test.{dataset}.out"
    log:
        "logs/a.{dataset}.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"