from typing import List

def parse_input(infile, parser, **kwargs):
    def inner(wildcards, input, output):
        with open(infile, "r") as fh:
            if parser is None:
                return fh.read().strip()
            else:
                return parser(fh, **kwargs)

    return inner


def extract_checksum(infile, **kwargs):
    try:
        import pandas as pd

        fix_file_name = lambda x: x.removeprefix("./")
        return (
            pd.read_csv(
                infile,
                sep="  ",
                header=None,
                engine="python",
                converters={1: fix_file_name},
            )
            .set_index(1)
            .loc[fix_file_name(kwargs.get("file"))]
            .item()
        )
    except ImportError:
        raise WorkflowError("Pandas is required to extract checksum from file.")
