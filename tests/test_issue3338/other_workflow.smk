rule all:
    input:
        "test_file.tsv"

rule test:
    output:
        file = "test_file.tsv"
    params:
        "pos1",
        "pos2",
        a = 3,
        b = 6
    shell:
        "echo {params} > {output.file}"

