from typing import List


def parse_input(infile, parser, **kwargs):
    def _parse_input(wildcards, input):
        _infile = infile(wildcards, input) if callable(infile) else infile
        with open(_infile, "r") as fh:
            if parser is None:
                return fh.read().strip()
            else:
                return parser(fh, **kwargs)

    return _parse_input


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


def flatten(list_of_lists: List) -> List:
    """Flatten an irregular list of lists recursively

    https://stackoverflow.com/a/53778278

    :param list_of_lists: A list of lists
    :return result: A list that has been flattened from a list of lists
    """
    result = list()
    for i in list_of_lists:
        if isinstance(i, list):
            result.extend(flatten(i))
        else:
            result.append(str(i))
    return result
