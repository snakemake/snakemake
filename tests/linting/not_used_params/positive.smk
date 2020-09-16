configfile: "config.yaml"

rule a:
    output:
        "{sample}.out"
    params:
        threshold=config["threshold"]
    log:
        "logs/{sample}.log"
    conda:
        "envs/software.yaml"
    shell:
        "echo {input} {wildcards.sample} {params.threshold} > {output} 2> {log}"


rule b:
    input:
        "test.in"
    output:
        "test.out"
    log:
        "logs/test.log"
    conda:
        "envs/software.yaml"
    shell:
        "awk {{print $1}} <{input} > {output} 2> {log}"
