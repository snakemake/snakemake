#!/bin/bash
set -e && echo "$0 $*" >&2

set -vx

sbatch_params=$1
dependencies=$2
script=$3

echo `bash slurm/sbatch.sh $sbatch_params $script "$dependencies"`
