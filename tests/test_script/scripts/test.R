print(snakemake@wildcards)
print(snakemake@threads)
print(snakemake@log)
print(snakemake@config)
print(snakemake@params)

if (!is.null(snakemake@params[["null_param"]])) {
    stop("Error evaluating null value.")
}
if (snakemake@params[["logical_param"]] != TRUE) {
    stop("Error evaluating logical value.")
}
if (snakemake@params[["integer_param"]] != 123L || typeof(snakemake@params[["integer_param"]]) != "integer") {
    stop("Error evaluating integer.")
}
if (snakemake@params[["double_param"]] != 123.0 || typeof(snakemake@params[["double_param"]]) != "double") {
    stop("Error evaluating double.")
}
if (!is.nan(snakemake@params[["nan_param"]])) {
    stop("Error evaluating NaN.")
}
if (!is.infinite(snakemake@params[["inf_param"]]) || snakemake@params[["inf_param"]] < 0) {
    stop("Error evaluating infinity.")
}
if (!is.infinite(snakemake@params[["neginf_param"]]) || snakemake@params[["neginf_param"]] > 0) {
    stop("Error evaluating negative infinity.")
}
if (snakemake@params[["complex_param"]] != 1+2i) {
    stop("Error evaluating complex.")
}
if (snakemake@params[["character_param"]] != "abc") {
    stop("Error evaluating character.")
}
if (!all.equal(snakemake@params[["vector_param"]], c(1, 2, 3))) {
    stop("Error evaluating vector.")
}
if (!identical(snakemake@params[["list_param"]], list(TRUE, 123L, "abc"))) {
    stop("Error evaluating list.")
}

if (snakemake@config[["test"]] != TRUE) {
    stop("Error evaluating config.")
}

if (snakemake@config[["foo\' bar"]] != "let\'s go") {
    stop("Error with the key/value containing single quotes.")
}

values <- scan(snakemake@input[[1]])
write(values, file = snakemake@output[["txt"]])
