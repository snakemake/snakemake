from pathlib import PosixPath

from snakemake.io import (
    WILDCARD_REGEX,
    AnnotatedString,
    IOFile,
    directory,
    expand,
    flag,
    is_flagged,
    temp,
)
from snakemake.exceptions import WildcardError, WorkflowError
from snakemake.path_modifier import PATH_MODIFIER_FLAG


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
    assert sorted(expand("{a}{b}{c}", **wildcards)) == sorted(
        ["135", "145", "235", "245"]
    )

    # redundant wildcards are provided
    assert sorted(expand("{a}{c}", **wildcards)) == sorted(["15", "25"])

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
        expand("first letter of sample: {samples[0]}", samples=["A123", "B456", "C789"])
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
        expand(
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
        expand(["a: {a} + b: {b}", "c: {c}"], a="aa", b=["b", "bb"], c=["c", "cc"])
    ) == sorted(["a: aa + b: b", "a: aa + b: bb", "c: c", "c: cc"])

    # expand on pathlib.Path objects
    assert expand(PosixPath() / "{x}" / "{y}", x="Hello", y="world") == ["Hello/world"]


def test_iofile_format_preserves_flags():
    """.format() on a flagged _IOFile (e.g. temp(), directory()) is an inherited plain
    str.format() call, which has no notion of flags and used to silently return a bare str
    with the flag gone entirely. It must instead carry the flags over to the formatted
    result, mirroring what apply_wildcards() already does for wildcard substitution."""
    flagged = IOFile(temp("results/{sample}.txt"), rule=None)

    formatted = flagged.format(sample="foo")
    assert formatted == "results/foo.txt"
    assert is_flagged(formatted, "temp")

    # a plain, unflagged file has nothing to preserve, so .format() keeps working as ordinary
    # str.format() -- this override must not affect the common, harmless case.
    unflagged = IOFile("results/{sample}.txt", rule=None)
    assert unflagged.format(sample="foo") == "results/foo.txt"
    assert not is_flagged(unflagged.format(sample="foo"), "temp")


def test_iofile_format_preserves_multiple_flags():
    flagged = IOFile(directory(temp("results/{sample}")), rule=None)

    formatted = flagged.format(sample="foo")
    assert formatted == "results/foo"
    assert is_flagged(formatted, "temp")
    assert is_flagged(formatted, "directory")


def test_iofile_format_preserves_path_modifier_flag():
    """Formatting must preserve PATH_MODIFIER_FLAG to avoid applying
    path modifiers twice."""
    flagged = IOFile(flag("results/{sample}.txt", PATH_MODIFIER_FLAG), rule=None)

    formatted = flagged.format(sample="foo")
    assert formatted == "results/foo.txt"
    assert is_flagged(formatted, PATH_MODIFIER_FLAG)


def test_iofile_format_raises_on_storage_file():
    """Storage flags cannot be copied after formatting because their
    query must remain consistent with the formatted local path."""

    class FakeStorageObject:
        def __init__(self, query):
            self.query = query

    flagged_storage = IOFile(
        flag(
            "results/{sample}.txt",
            "storage_object",
            FakeStorageObject("s3://bucket/{sample}.txt"),
        ),
        rule=None,
    )

    try:
        flagged_storage.format(sample="foo")
        raise AssertionError(
            ".format() on a storage-flagged _IOFile should raise WorkflowError"
        )
    except WorkflowError as e:
        # _IOFile has an apply_wildcards() alternative that AnnotatedString lacks, so
        # it must get the more specific message pointing users to it.
        assert "apply_wildcards" in str(e)


def test_annotated_string_format_preserves_flags():
    """temp(), directory(), etc. attach flags to a plain AnnotatedString before it is
    ever wrapped as an _IOFile, so .format() can be called on it directly."""
    flagged = temp("results/{sample}.txt")
    assert isinstance(flagged, AnnotatedString)

    formatted = flagged.format(sample="foo")
    assert formatted == "results/foo.txt"
    assert is_flagged(formatted, "temp")

    # a plain, unflagged AnnotatedString has nothing to preserve, so .format() keeps
    # working as ordinary str.format() -- this override must not affect that case.
    unflagged = AnnotatedString("results/{sample}.txt")
    assert unflagged.format(sample="foo") == "results/foo.txt"
    assert not is_flagged(unflagged.format(sample="foo"), "temp")


def test_annotated_string_format_preserves_multiple_flags():
    flagged = directory(temp("results/{sample}"))
    assert isinstance(flagged, AnnotatedString)

    formatted = flagged.format(sample="foo")
    assert formatted == "results/foo"
    assert is_flagged(formatted, "temp")
    assert is_flagged(formatted, "directory")


def test_annotated_string_format_raises_on_storage_flag():
    """Mirrors test_iofile_format_raises_on_storage_file(), but for the AnnotatedString
    that storage() returns directly, before it is ever wrapped in an _IOFile."""

    class FakeStorageObject:
        def __init__(self, query):
            self.query = query

    flagged_storage = flag(
        "results/{sample}.txt",
        "storage_object",
        FakeStorageObject("s3://bucket/{sample}.txt"),
    )
    assert isinstance(flagged_storage, AnnotatedString)

    try:
        flagged_storage.format(sample="foo")
        raise AssertionError(
            ".format() on a storage-flagged AnnotatedString should raise WorkflowError"
        )
    except WorkflowError as e:
        # AnnotatedString has no apply_wildcards(), so it must not get the _IOFile-
        # specific message that points users to it.
        assert "apply_wildcards" not in str(e)
