#!python

import argparse
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument("--name")
parser.add_argument("--format")
parser.add_argument("-X", action="store_true")
parser.add_argument("--parsable2", action="store_true")
parser.add_argument("--noheader", action="store_true")

args = parser.parse_args()

sp.call(["squeue", "--noheader", "--format", "%F|%T", "--name", args.name])