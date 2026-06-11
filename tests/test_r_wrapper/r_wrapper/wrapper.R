outfile = snakemake@output[[1]]
x <- data.frame()
write.table(x, file=outfile, col.names=FALSE)
