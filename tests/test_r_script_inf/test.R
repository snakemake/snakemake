#! /opt/conda/envs/dev/bin/R

run <- function (output, threshold){
  log <- c(
    paste("The type: ", class(threshold[1]), ".", sep = ""),
    paste("The value: ", paste(threshold, collapse = ', '), ".", sep = "")
  )
  conn <- file(output, "w")
  writeLines(log, conn)
  close(conn)
}

run(
  snakemake@output[[1]],
  snakemake@params[["threshold"]]
)