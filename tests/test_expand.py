from snakemake.io import expand
from snakemake.exceptions import WildcardError
import pytest


def test_simple_expand():
    # single filepattern
    assert expand("{a}.out", a="test") == ["test.out"]
    # multiple filepatterns
    assert expand(["{a}.out", "{b}.out"], a="a", b="b") == ["a.out", "b.out"]
    # multiple wildcards
    assert expand("{a}.out", a=["1", "2", "3"]) == ["1.out", "2.out", "3.out"]
    # multiple wildcards and patterns
    assert expand(["{a}_{b}.ab", "{b}.b"], a="1 2".split(), b="3 4".split()) == [
        "1_3.ab",
        "1_4.ab",
        "2_3.ab",
        "2_4.ab",
        "3.b",
        "4.b",
    ]
    # replace product
    assert expand(["{a}_{b}.ab", "{b}.b"], zip, a="1 2".split(), b="3 4".split()) == [
        "1_3.ab",
        "2_4.ab",
        "3.b",
        "4.b",
    ]


def test_allow_missing():
    # single filepattern
    assert expand("{a}_{b}.out", allow_missing=True) == ["{a}_{b}.out"]
    assert expand("{a}_{b}.out", a="test", allow_missing=True) == ["test_{b}.out"]
    # none missing
    assert expand("{a}.out", a="test", allow_missing=True) == ["test.out"]
    # wildcard is allow_missing
    assert expand("{allow_missing}.out", allow_missing=True) == ["True.out"]
    # allow_missing not True
    assert expand("{a}.out", a="test", allow_missing="test2") == ["test.out"]
    with pytest.raises(WildcardError) as e:
        expand("{a}.out", allow_missing="test2")
    assert str(e.value) == "No values given for wildcard 'a'."

    # multiple filepatterns
    assert expand(["{a}.out", "{b}.out"], allow_missing=True) == ["{a}.out", "{b}.out"]
    # multiple wildcards
    assert expand("{a}_{b}.out", a=["1", "2", "3"], allow_missing=True) == [
        "1_{b}.out",
        "2_{b}.out",
        "3_{b}.out",
    ]
    # multiple wildcards and patterns
    assert expand(
        ["{a}_{b}_{C}.ab", "{b}_{c}.b"],
        a="1 2".split(),
        b="3 4".split(),
        allow_missing=True,
    ) == ["1_3_{C}.ab", "1_4_{C}.ab", "2_3_{C}.ab", "2_4_{C}.ab", "3_{c}.b", "4_{c}.b"]
    # replace product
    assert expand(
        ["{a}_{b}_{C}.ab", "{b}_{c}.b"],
        zip,
        a="1 2".split(),
        b="3 4".split(),
        allow_missing=True,
    ) == ["1_3_{C}.ab", "2_4_{C}.ab", "3_{c}.b", "4_{c}.b"]
