import snakemake.io


# Alias for expand that provides a more intuitive name for the use case of
# collecting files from previous jobs.
collect = snakemake.io.expand
