rule a:
    input:
        "test.in"
    output:
        "test.out"
    log:
        "logs/a.log"
    singularity:
        "docker://bash"
    shell:
        "cp {input} {output}"