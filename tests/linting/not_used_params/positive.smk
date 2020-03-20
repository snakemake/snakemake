configfile: "config.yaml"

rule a:
    output:
        "test.out"
    params:
        threshold=config["threshold"]
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "echo {params.threshold} > {output}"