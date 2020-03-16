rule a:
    input:
        "test.in"
    output:
        "test.out"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"