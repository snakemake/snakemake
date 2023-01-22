#!/usr/bin/env bash

echo "The first input file is ${snakemake_input[0]}" > "${snakemake_output[0]}" 2> "${snakemake_log[0]}"
echo "The named input file is ${snakemake_input[named]}" >> "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"
echo "The requested number of threads is ${snakemake[threads]}" >> "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"

