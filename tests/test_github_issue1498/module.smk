rule a:
    output:
        a=touch("results/a.txt"),


def unpack_a(wildcards):
    return {"a": rules.a.output.a}


def lambda_a(wildcards):
    return rules.a.output.a


rule b:
    input:
        unpack(unpack_a),
        # a = lambda_a,
    output:
        b="results/b.txt",
    shell:
        "cp {input.a} {output}"
