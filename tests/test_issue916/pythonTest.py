#!/usr/bin/env python
import sys
import os
import local_script

# Ensure that the path ends in .snakemake
dname = os.path.dirname(__file__)
print(dname)
if not dname.endswith(".snakemake/scripts"):
    sys.exit("We're not being written in the output directory!\n")
of = open(snakemake.output[0], "w")
of.write("local_script.contents {}\n".format(local_script.contents))
of.close()
