rule all:
    input:
        expand("results/{i}.txt", i=range(5)),


rule a:
    output:
        "results/{i}.txt",
    shell:
        "echo {wildcards.i} > {output}"
