#!/usr/bin/env python3
import sys
import os
import random

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)
with open("qsub.log", "a") as log:
    print(job_properties, file=log)

print(random.randint(1, 100))
os.system("sh {}".format(jobscript))
