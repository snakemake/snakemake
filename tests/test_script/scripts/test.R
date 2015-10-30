print(snakemake@wildcards)
print(snakemake@threads)
print(snakemake@log)
print(snakemake@config)
print(snakemake@params)

values <- scan(snakemake@input[[1]])
write(values, file = snakemake@output[["txt"]])
