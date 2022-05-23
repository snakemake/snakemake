import pandas as pd
import pytest
from snakemake.utils import Paramspace

MOCKDICT = {
    "first": [1, 2, 3],
    "second": "A B C".split(),
    "third.a": "A B C".split(),
    "third.c": [4, 5, 6],
    "four": "AA BA CA".split(),
    "five_or_six.2": 3 * ["a"],
}


@pytest.fixture
def mockdf():
    df = pd.DataFrame.from_dict(MOCKDICT)
    return df


def test_paramspace_constructs(mockdf):

    spc = Paramspace(mockdf)
    assert len(spc.dataframe) == 3


def test_instances_have_keys(mockdf):

    spc = Paramspace(mockdf)
    obs = list(spc.instance_patterns)
    assert len(obs) > 0
    assert len(obs) == len(mockdf)

    for coln in mockdf.columns:
        assert coln in obs[0]
