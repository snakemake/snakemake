rule a:
    input:
        "test.in"
    output:
        "test.out"
    log:
        "logs/a.log"
    shell:
        "cp {input} {output}"