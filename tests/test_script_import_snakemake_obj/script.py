from snakemake.script import snakemake, is_script

if not is_script:
    raise Exception("This script can only be run as a Snakemake script.")

with open(snakemake.output[0], "w") as f:
    f.write("Hello world!")
