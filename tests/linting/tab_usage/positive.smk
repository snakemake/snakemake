rule a:
    input:
        "test.txt"
    output:
        "test.out"
    log:
        "logs/a.log"
    conda:
        "envs/a.yaml"
    shell:
        "cat {input} > {output}"
