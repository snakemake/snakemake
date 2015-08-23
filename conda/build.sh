#!/bin/sh

conda build snakemake
conda convert ~/miniconda3/conda-bld/linux-64/snakemake-*.tar.bz2 -p all
echo Due to a bug in conda convert, the subdir in info/index.json is not updated (to e.g. win-64).
echo This has to be done manually in the tarball.
exit 0
# TODO reactivate this once conda convert has been fixed.
for p in */snakemake-*.tar.bz2
do
    binstar upload $p
done
