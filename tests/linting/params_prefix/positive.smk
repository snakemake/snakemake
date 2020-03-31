rule a:
    input:
        "path/test.txt"
    output:
        "test.txt"
    params:
        prefix=lambda w, input: os.path.dirname(input[0])
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"