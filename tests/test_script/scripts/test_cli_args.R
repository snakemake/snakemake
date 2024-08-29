
library(docopt)
'test_cli_args.

Usage:
  test_cli_args.R <message>
' -> doc


arguments <- docopt(doc)
print(arguments)


write(arguments$message, file = snakemake@output[["txt"]])
