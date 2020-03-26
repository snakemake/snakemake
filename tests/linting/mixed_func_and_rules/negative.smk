def get_dataset(wildcards):
    return config["datasets"]

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