#!/bin/sh

conda build snakemake
conda convert ~/miniconda3/conda-bld/linux-64/snakemake-*.tar.bz2 -p all
binstar upload --force */snakemake-*.tar.bz2
