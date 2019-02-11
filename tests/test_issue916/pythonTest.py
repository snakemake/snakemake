#!/usr/bin/env python
import sys
import os
import local_script
from pathlib import Path

# Ensure that the __real_file__ path ends in .snakemake
dname = Path(__real_file__).parent
print(dname)
if not dname.parts[-2:] == (".snakemake", "scripts"):
    sys.exit("We're not being written in the output directory!\n")

# Ensure that the __file__ path ends in test_issue916
dname = Path(__file__).parent
print(dname)
if not dname.name == "test_issue916":
    sys.exit("We're not faking __file__ properly!\n")

# Write out script to indicte success.
of = open(snakemake.output[0], "w")
of.write("local_script.contents {}\n".format(local_script.contents))
of.close()
