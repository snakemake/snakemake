from snakemake import shell
shell.use_bash_on_win()

shell("echo {} > {}".format(snakemake.params['title'], snakemake.output[0]))
