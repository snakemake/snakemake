#!/bin/bash
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

# setup conda
if type conda > /dev/null; then exit 0; fi
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p miniconda
conda env create --quiet --name snakemake --file test-environment.yml
conda create --quiet -n black black

