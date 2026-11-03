rule produce:
    output:
        "test.out"
    shell:
        "echo test > {output}"
