#!/bin/bash
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
source .circleci/common.sh

sudo add-apt-repository -y ppa:gophers/archive
sudo apt-get update

# install stress
sudo apt-get install stress

# install singularity
sudo apt-get install build-essential libssl-dev uuid-dev libgpgme11-dev libseccomp-dev wget pkg-config squashfs-tools libarchive-dev golang-1.11
export PATH=/usr/lib/go-1.11/bin:$PATH
wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VER}/singularity-${SINGULARITY_VER}.tar.gz
tar -xvf singularity-$SINGULARITY_VER.tar.gz
cd singularity
./mconfig
make -C builddir
sudo make -C builddir install
