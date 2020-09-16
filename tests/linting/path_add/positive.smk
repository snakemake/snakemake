rule a:
    input:
        "test.in"
    output:
        "results/tables/test.out"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"