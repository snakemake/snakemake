print(snakemake@wildcards)
print(snakemake@threads)
print(snakemake@log)
print(snakemake@config)
print(snakemake@params)
if(snakemake@params[["xy"]] != TRUE) {
    stop("Error evaluating param.")
}
if(snakemake@config[["test"]] != TRUE) {
    stop("Error evaluating config.")
}

values <- scan(snakemake@input[[1]])
write(values, file = snakemake@output[["txt"]])
