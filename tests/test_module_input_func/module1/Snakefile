def txt_output(wildcards):
    return ["results/C.txt"]


rule all:
    input:
        txt_output


rule txt:
    output:
        "results/C.txt"
    shell:
        "echo 'C' "
        ">{output} "
