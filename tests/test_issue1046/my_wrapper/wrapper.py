from snakemake import shell

shell("echo {} > {}".format(snakemake.params["title"], snakemake.output[0]))
