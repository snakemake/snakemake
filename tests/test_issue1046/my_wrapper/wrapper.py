from snakemake import shell
shell.executable("bash")

shell("echo {} > {}".format(snakemake.params["title"], snakemake.output[0]))
