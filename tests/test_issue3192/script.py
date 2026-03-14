import shutil

shutil.copyfile(
    snakemake.input[0],
    snakemake.output[0])
