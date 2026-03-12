#!/usr/bin/env python3

try:
    from .another_module import another_function
    from ..helpers import helper_function
except ImportError:
    print("Relative imports do not work correctly. Aborting.")
    exit()


def main():
    print("Relative imports work correctly")
    another_function()
    helper_function()


if __name__ == "__main__":
    main()
