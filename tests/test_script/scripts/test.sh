#!/usr/bin/env bash

echo "The first input file is ${SNAKEMAKE_INPUT[0]}" > "${SNAKEMAKE_OUTPUT[0]}" 2> "${SNAKEMAKE_LOG[0]}"
echo "The named input file is ${SNAKEMAKE_INPUT[named]}" >> "${SNAKEMAKE_OUTPUT[0]}" 2>> "${SNAKEMAKE_LOG[0]}"
echo "The requested number of threads is ${SNAKEMAKE[threads]}" >> "${SNAKEMAKE_OUTPUT[0]}" 2>> "${SNAKEMAKE_LOG[0]}"

