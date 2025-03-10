rule all:
    input:
        "test_1_file.tsv",
        "test_2_file.tsv"

rule test:
    input:
        file = "file_{nr}.txt"
    output:
        file = "test_{nr}_file.tsv"
    params:
        a = 3,
        b = 6
    shell:
        "echo {params.a} > {output.file} && echo {params.b} >> {output.file}"

