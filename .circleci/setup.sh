#!/bin/bash
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

# conda
if type conda > /dev/null; then exit 0; fi
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p miniconda
conda env create --name snakemake --file test-environment.yml
conda create -n black black

# stress
sudo apt-get update
sudo apt-get install stress

# singularity
source .circleci/common.sh
sudo add-apt-repository -y ppa:gophers/archive
sudo apt-get update
sudo apt-get install build-essential libssl-dev uuid-dev libgpgme11-dev libseccomp-dev wget pkg-config squashfs-tools libarchive-dev golang-1.11
export PATH=/usr/lib/go-1.11/bin:$PATH
wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VER}/singularity-${SINGULARITY_VER}.tar.gz
tar -xvf singularity-$SINGULARITY_VER.tar.gz
cd singularity
./mconfig
make -C builddir
sudo make -C builddir install
