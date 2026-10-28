#!/usr/bin/env python3

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--input", action="append")

args = parser.parse_args()

print("Input: ", args.input)
