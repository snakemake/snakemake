os.environ["THRESHOLD"] = "1.0"

envvars:
    "THRESHOLD"

rule a:
    output:
        "test.out"
    params:
        threshold=os.environ["THRESHOLD"]
    log:
        "logs/a.log"
    conda:
        "envs/software.yaml"
    shell:
        "echo {params.threshold} > {output}"