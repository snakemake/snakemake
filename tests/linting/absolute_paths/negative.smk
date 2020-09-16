rule a:
    input:
        "/home/johannes/someproject/test.in"
    output:
        "test.out"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"