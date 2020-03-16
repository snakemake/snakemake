rule a:
    input:
        "test.in"
    output:
        dynamic("test.out")
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"