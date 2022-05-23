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


def test_instance_patterns_have_keys(mockdf):

    spc = Paramspace(mockdf)
    obs = list(spc.instance_patterns)
    assert len(obs) > 0
    assert len(obs) == len(mockdf)

    for coln in mockdf.columns:
        assert coln in obs[0]


def test_instance_patterns_contain_correct_rows(mockdf):

    spc = Paramspace(mockdf)
    nrows = len(mockdf)
    ncols = len(mockdf.columns)

    observed = list(spc.instance_patterns)
    colnames = list(spc.ordered_columns)

    for rowi in range(nrows):

        obs = observed[rowi]
        entries = obs.split("/")
        assert len(entries) == ncols

        for coli in range(ncols):
            pair = entries[coli]
            colname = colnames[coli]

            assert colname in pair
            assert str(mockdf[colname][rowi]) in pair
