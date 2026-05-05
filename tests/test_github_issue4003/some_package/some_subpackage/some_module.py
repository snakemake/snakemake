#!/usr/bin/env python3

import contextlib
import sys


try:
    output = open(snakemake.output["stdout"], "w")
except NameError:
    output = sys.stdout

try:
    from .another_module import another_function
    from ..helpers import helper_function
except ImportError:
    print("Relative imports do not work correctly. Aborting.", file=output)
    output.close()
    exit()


def main():
    with contextlib.redirect_stdout(output):
        print("Relative imports work correctly")
        another_function()
        helper_function()

    output.close()


if __name__ == "__main__":
    main()
