rule a:
    input:
        "test.in"
    output:
        "test.out"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "cp {input} {output}"


rule b:
    input:
        "test.in"
    output:
        "test2.out"
    log:
        "logs/a.log"
    container:
        "docker://alpine"
    shell:
        "cp {input} {output}"