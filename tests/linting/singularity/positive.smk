rule a:
    input:
        "test.in"
    output:
        "test.out"
    log:
        "logs/a.log"
    container:
        "docker://bash"
    shell:
        "cp {input} {output}"