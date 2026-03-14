#!/usr/bin/env bash
set -euo pipefail
exec > "${snakemake_output[0]}" 2> "${snakemake_log[0]}"
set -x

# Awkward characters are robustly quoted?
test "${snakemake_params[astring]}" = $'foo\n\'\\\" ' || echo MISMATCH

echo "The first input file is ${snakemake_input[0]}"
echo "The named input file is ${snakemake_input[named]}"
echo "The requested number of threads is ${snakemake[threads]}"

echo "snakemake_config is type ${snakemake_config@a}"

echo "The list passed as a parameter is *${snakemake_params[alist]}*"
echo "The wildcards are *${snakemake_wildcards[bash]}* and *${snakemake_wildcards[empty]}*"
echo "The config items are *${snakemake_config[test]}* *${snakemake_config[testint]}* *${snakemake_config[testfloat]}*"
echo "The config item with quotes in is *${snakemake_config['foo'"'"' bar']}*"
