if not snakemake.resources.mem_mb == 20:
    raise ValueError("Handover of all resources did not work.")
with open(snakemake.output[0], "w") as out:
    print("test", file=out)