#!/usr/bin/env python3
import sys
import os
import random
import re

jobscript = sys.argv[1]

os.system("strace -f -o trace {}".format(jobscript))

with open("trace", "r") as fp:
    lines = fp.readlines()

regex = re.compile('.*openat\(.*.snakemake/incomplete".*')
matches = list(filter(lambda l: regex.match(l), lines))
if len(matches) > 0:
    sys.stderr.write(repr(matches))
    sys.exit(1)

print(random.randint(1, 100))
