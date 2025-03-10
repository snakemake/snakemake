import json

with open(snakemake.output[0], "w") as out:
    # convert numpy types to native python types for serialization
    json.dump({k: v.item() for k, v in snakemake.params.simulation.items()}, out)