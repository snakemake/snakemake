import json

with open(snakemake.output[0], "w") as out:
    json.dump(snakemake.params.simulation, out)