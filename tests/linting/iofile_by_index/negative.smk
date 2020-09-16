rule a:
    input:
        "test.in",
        "test.aux"
    output:
        "test.out"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input[0]} {output}"