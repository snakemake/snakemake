echo @($(bcftools --help).splitlines()[1]) > @(snakemake.output[0])
