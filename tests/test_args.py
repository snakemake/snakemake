__authors__ = ["K.D. Murray"]
__copyright__ = "Copyright 2024, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


import pytest

from snakemake.settings.types import Batch


def test_parse_batch():
    from snakemake.cli import parse_batch

    assert parse_batch("aggregate=1/2") == Batch("aggregate", 1, 2)


def test_parse_config():
    from snakemake.cli import parse_config

    assert parse_config(["k=24", "w=1000", "_=x", "ab=1", "key-with-hyphens=a"]) == {
        "k": 24,
        "w": 1000,
        "_": "x",
        "ab": 1,
        "key-with-hyphens": "a",
    }
    for key in ["1k", "-k", "k-", "-"]:
        with pytest.raises(ValueError):
            parse_config([f"{key}=1"])
