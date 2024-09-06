from snakemake.script import snakemake

with open(snakemake.output[0], "w") as out:
    print(1, 2, 3, file=out)
