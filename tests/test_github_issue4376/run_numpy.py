import numpy

with open(snakemake.output[0], "w") as f:
    f.write(numpy.__version__)
