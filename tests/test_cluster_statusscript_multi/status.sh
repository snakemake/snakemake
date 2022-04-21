#!/bin/bash

# The argument passed from sbatch is "jobid;cluster_name"

arg="$1"
jobid="${arg%%;*}"
cluster="${arg##*;}"

echo success
