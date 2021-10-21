
infile <- snakemake@input[["text_file"]]

# Check the read value is as expected in R.
param_dir <- snakemake@params[["param_dir"]]
stopifnot(param_dir == "dir")

outfile <- snakemake@output[["out_file"]]
file.create(outfile)
