from snakemake.shell import shell

shell("echo Hello World > {snakemake.output}")
