rule a:
    input:
        "test.in"
    output:
        "test.out"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    version: "1.0"
    shell:
        "cp {input} {output}"