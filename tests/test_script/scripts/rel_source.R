snakemake@source("source_me.R")
cat(hi(),file=snakemake@output[[1]])
