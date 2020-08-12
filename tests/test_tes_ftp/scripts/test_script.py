with open(snakemake.output[0], "w") as f:
  f.write(snakemake.params.message)
