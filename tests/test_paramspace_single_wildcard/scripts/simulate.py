with open(snakemake.output[0], "w") as out:
    print(snakemake.params.simulation, file=out)