rule a:
    input:
        some="test.in",
        other="test.aux"
    output:
        "test.out"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input.some} {output}"