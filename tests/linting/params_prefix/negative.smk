rule a:
    input:
        "path/test.txt"
    output:
        "test.txt"
    params:
        prefix="path"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"