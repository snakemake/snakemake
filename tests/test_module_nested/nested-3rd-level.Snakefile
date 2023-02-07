from snakemake.utils import min_version
min_version("6.0")

rule leaf:
    output:
        ".done"
    shell:
        "touch {output}"

rule default:
    input:
        rules.leaf.output
    default_target: True
