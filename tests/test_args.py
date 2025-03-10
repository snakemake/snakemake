__authors__ = ["K.D. Murray"]
__copyright__ = "Copyright 2024, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


from snakemake.settings.types import Batch


def test_parse_batch():
    from snakemake.cli import parse_batch

    assert parse_batch("aggregate=1/2") == Batch("aggregate", 1, 2)
