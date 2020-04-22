#!/usr/bin/env python3
import sys
import os
import random

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

expected_model = "nvidia-tesla-p100" if job_properties["rule"] == "a" else "nvidia-tesla-1000"
print(job_properties)

assert job_properties["resources"]["gpu_model"] == expected_model
with open("qsub.log", "a") as log:
    print(job_properties, file=log)

print(random.randint(1, 100))
os.system("sh {}".format(jobscript))
