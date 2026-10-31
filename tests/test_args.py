__authors__ = ["K.D. Murray"]
__copyright__ = "Copyright 2024, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


import pytest

from snakemake.exceptions import WorkflowError
from snakemake.settings.types import Batch, MaxJobsPerTimespan


def test_parse_batch():
    from snakemake.cli import parse_batch

    assert parse_batch("aggregate=1/2") == Batch("aggregate", 1, 2)


@pytest.mark.parametrize(
    "arg,max_jobs,timespan",
    [
        ("100/1s", 100, 1),
        ("50/1m", 50, 60),
        ("10/2h", 10, 7200),
        ("100/1ms", 100, 0.001),
        ("0.5/1s", 1, 2),
        ("2.5/1m", 5, 120),
    ],
)
def test_parse_max_jobs_per_timespan(arg, max_jobs, timespan):
    limit = MaxJobsPerTimespan.parse_choice(arg)
    assert limit.max_jobs == max_jobs
    assert isinstance(limit.max_jobs, int)
    assert limit.timespan == pytest.approx(timespan)


@pytest.mark.parametrize("arg", ["1/1x", "1/1msx", "1/s", "a/1s", ".5/1s", "1/1s "])
def test_parse_max_jobs_per_timespan_invalid(arg):
    with pytest.raises(WorkflowError):
        MaxJobsPerTimespan.parse_choice(arg)
