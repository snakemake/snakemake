from pathlib import PosixPath

from snakemake.io import WILDCARD_REGEX, expand
from snakemake.exceptions import WildcardError


def test_wildcard_regex():
    def matches(text):
        return [
            (match.group("name"), match.group("constraint"))
            for match in WILDCARD_REGEX.finditer(text)
        ]

    # without constraints
    assert matches("") == []
    assert matches("{") == []
    assert matches("}") == []
    assert matches("{}") == []
    assert matches("{0}") == [("0", None)]
    assert matches("{abc}") == [("abc", None)]
    assert matches("abc{def}{ghi}") == [("def", None), ("ghi", None)]

    # with constraints
    assert matches("{w,constraint}") == [("w", "constraint")]
    assert matches("{w , constraint}") == [("w", "constraint")]
    # fails because constraint is detected as 'constraint '
    # assert matches('{w,constraint }') == [('w', 'constraint')]
    assert matches("abc { w , constraint} def") == [("w", "constraint")]

    # multiple wildcards
    assert matches("{a,1} {b,2} {c,3}") == [("a", "1"), ("b", "2"), ("c", "3")]

    # more complicated constraints
    assert matches(r"{w,([a-z]+|pat\|t*ern)}") == [("w", r"([a-z]+|pat\|t*ern)")]
    assert matches(r"{w,([a-z]+|pat\|te{1,3}rn){5,7}}") == [
        ("w", r"([a-z]+|pat\|te{1,3}rn){5,7}")
    ]

    # This used to be very slow with an older version of the regex
    assert matches("{w, long constraint without closing brace") == []


def test_expand():
    wildcards = {"a": [1, 2], "b": [3, 4], "c": [5]}

    # each provided wildcard is used in the filepattern
    assert sorted(expand("{a}{b}{c}", **wildcards)) == sorted(  # type: ignore[reportArgumentType]
        ["135", "145", "235", "245"]
    )

    # redundant wildcards are provided
    assert sorted(expand("{a}{c}", **wildcards)) == sorted(["15", "25"])  # type: ignore[reportArgumentType]

    # missing wildcards (should fail)
    try:
        expand("{a}{d}", **wildcards)
        assert False
    except WildcardError:
        pass

    # do not expand on strings and non iterables
    assert expand("{x}{y}", **{"x": "Hello, ", "y": "world!"}) == ["Hello, world!"]
    assert expand("{x}{y}", **{"x": 4, "y": 2}) == ["42"]

    # format-minilang: field names
    assert sorted(
        expand("first letter of sample: {samples[0]}", samples=["A123", "B456", "C789"])  # type: ignore[reportArgumentType]
    ) == sorted(
        [
            "first letter of sample: A",
            "first letter of sample: B",
            "first letter of sample: C",
        ]
    )
    assert expand("{str.__class__}", str="") == ["<class 'str'>"]

    # format-minilang: conversions
    class ConvTest:
        def __str__(self):
            return "string"

        def __repr__(self):
            return "representation"

    assert expand("{test!r}", test=ConvTest()) == ["representation"]
    assert expand("{test!s}", test=ConvTest()) == ["string"]

    # format-minilang: format specifications
    assert sorted(
        expand(  # type: ignore[reportArgumentType]
            "The answer to life, the universe, and everything: {answer:f}",
            answer=range(41, 43),
        )
    ) == sorted(
        [
            "The answer to life, the universe, and everything: 41.000000",
            "The answer to life, the universe, and everything: 42.000000",
        ]
    )

    # multiple filepatterns with different wildcards
    assert sorted(
        expand(["a: {a} + b: {b}", "c: {c}"], a="aa", b=["b", "bb"], c=["c", "cc"])  # type: ignore[reportArgumentType]
    ) == sorted(["a: aa + b: b", "a: aa + b: bb", "c: c", "c: cc"])

    # expand on pathlib.Path objects
    assert expand(PosixPath() / "{x}" / "{y}", x="Hello", y="world") == ["Hello/world"]


def test_typed():
    from snakemake.io.typed import typed_factory
    from tempfile import TemporaryDirectory
    from typing import NamedTuple, List

    class CustomType(NamedTuple):
        a: int
        b: List[int]

    with TemporaryDirectory() as tmpdir:
        typed = typed_factory(CustomType, loader="json")
        file = f"{tmpdir}/test.json"
        objd = {"a": 1, "b": [1, 2, 3]}
        obj = CustomType(**objd)  # type: ignore[arg-type]
        typed(file).dump(**objd)
        assert typed(file).load() == obj
        typed = typed_factory(CustomType, loader="yaml")
        file = f"{tmpdir}/test.yaml"
        typed(file).dump(**objd)
        assert typed(file).load() == obj
        typed = typed_factory(CustomType, loader="toml")
        file = f"{tmpdir}/test.toml"
        typed(file).dump(**objd)
        assert typed(file).load() == obj
        typed = typed_factory("pkl")
        file = f"{tmpdir}/test.pkl"
        typed(file).dump(objd)
        assert typed(file).load() == objd
        typed = typed_factory("npy")
        file = f"{tmpdir}/test.npy"
        import numpy as np

        arr = np.array([1, 2, 3])
        typed(file).dump(arr)
        assert np.array_equal(typed(file).load(), arr)
        import pandas as pd

        file = f"{tmpdir}/test.csv"
        typed = typed_factory(loader=pd.read_csv, dumper=pd.DataFrame.to_csv)  # type: ignore
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        typed(file).dump(df)
        assert typed(file).load().drop(columns="Unnamed: 0").equals(df)  # type: ignore[reportAttributeAccessIssue]
        typed = typed_factory("tsv")
        file = f"{tmpdir}/test.tsv"
        typed(file).dump(df)
        assert typed(file).load().equals(df)  # type: ignore[reportAttributeAccessIssue]
