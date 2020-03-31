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