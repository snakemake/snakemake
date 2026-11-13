with open(snakemake.output[0], "w") as out:
    print("script from local git module", file=out)
