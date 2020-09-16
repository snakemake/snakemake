configfile: "config.yaml"

rule a:
    output:
        "test.out"
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "echo {config[threshold]} > {output}"