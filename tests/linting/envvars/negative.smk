os.environ["THRESHOLD"] = "1.0"

rule a:
    input:
        "test.in"
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