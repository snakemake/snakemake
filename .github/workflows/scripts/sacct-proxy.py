#!python

import argparse
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument("--name")
parser.add_argument("--format")
parser.add_argument("-X", action="store_true")
parser.add_argument("--parsable2", action="store_true")
parser.add_argument("--noheader", action="store_true")
parser.add_argument("-n", action="store_true")
parser.add_argument("-u")
parser.add_argument("-o")


args = parser.parse_args()

if args.n and args.u and args.o:
    # mimic account query from the executor
    print("runner")
elif args.name:
    sp.call(["squeue", "--noheader", "--format", "%F|%T", "--name", args.name])
else:
    raise ValueError("Unsupported arguments")
