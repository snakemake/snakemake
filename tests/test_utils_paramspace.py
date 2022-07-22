from pathlib import Path

import pandas as pd
import pytest
from snakemake.io import expand
from snakemake.utils import Paramspace

MOCKDICT = {
    "first": [1, 2, 3],
    "second": "A B C".split(),
    "third.a": "A B C".split(),
    "third.c": [4, 5, 6],
    "four": "AA BA CA".split(),
    "five_or_six.2": 3 * ["a"],
}

PARAMS_TSV = {"alpha": [1, 2], "beta": [0.1, 0.0], "gamma": [0.99, 3.9]}
SPECIAL_TSV = {
    "alpha": [1, 2],
    "beta": [0.1, 0.0],
    "gamma": [0.99, 3.9],
    "delta": ["foo", "Bar"],
}


@pytest.fixture
def mockdf():
    df = pd.DataFrame.from_dict(MOCKDICT)
    return df


@pytest.fixture
def df():
    df_ = pd.DataFrame.from_dict(PARAMS_TSV)
    return df_


@pytest.fixture
def sdf():
    df_ = pd.DataFrame.from_dict(SPECIAL_TSV)
    return df_


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


def test_paramspace_one_simulation_beta(df):

    spc = Paramspace(df, filename_params=["beta"])
    obs = expand(
        "results/one/simulations/{params}.tsv",
        params=spc.instance_patterns,
    )

    assert obs[0] == "results/one/simulations/alpha~1/gamma~0.99/beta~0.1.tsv"


def test_paramspace_one_simulation_beta_nonstrict(df):

    spc = Paramspace(df, filename_params=["beta"], types_strict=False)
    obs = expand(
        "results/one/simulations/{params}.tsv",
        params=spc.instance_patterns,
    )

    assert obs[0] == "results/one/simulations/alpha~1.0/gamma~0.99/beta~0.1.tsv"


def test_paramspace_one_simulation_beta_wcpattern(df):

    spc = Paramspace(df, filename_params=["beta"])

    assert spc.wildcard_pattern != None
    assert spc.instance != None


def test_paramspace_one_simulation_delta_special(sdf):

    spc = Paramspace(sdf, filename_params=["delta"])
    obs = expand(
        "results/one/simulations/{params}.tsv",
        params=spc.instance_patterns,
    )

    assert obs[0] == "results/one/simulations/alpha~1/beta~0.1/gamma~0.99/delta~foo.tsv"


def test_paramspace_one_simulation_delta_special_nonstrict(sdf):

    spc = Paramspace(sdf, filename_params=["delta"], types_strict=False)
    obs = expand(
        "results/one/simulations/{params}.tsv",
        params=spc.instance_patterns,
    )

    assert obs[0] == "results/one/simulations/alpha~1/beta~0.1/gamma~0.99/delta~foo.tsv"


def test_paramspace_two_simulation_beta_gamma(df):

    spc = Paramspace(df, filename_params=["beta", "gamma"])
    obs = expand(
        "results/two/simulations/{params}.tsv",
        params=spc.instance_patterns,
    )

    assert obs[0] == "results/two/simulations/alpha~1/beta~0.1_gamma~0.99.tsv"
    assert obs[1] == "results/two/simulations/alpha~2/beta~0.0_gamma~3.9.tsv"


def test_paramspace_full_reorder_simulation_beta(df):

    spc = Paramspace(df, filename_params=["beta", "gamma", "alpha"])

    obs = expand(
        "results/full_reorder/simulations/{params}.tsv",
        params=spc.instance_patterns,
    )

    assert obs[0] == "results/full_reorder/simulations/beta~0.1_gamma~0.99_alpha~1.tsv"
    assert obs[1] == "results/full_reorder/simulations/beta~0.0_gamma~3.9_alpha~2.tsv"
