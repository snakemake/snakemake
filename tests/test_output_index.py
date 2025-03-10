from unittest.mock import Mock

import pytest

from snakemake.output_index import OutputIndex


@pytest.fixture()
def rule_0(mocker):
    mock_products_rule_0 = [
        mocker.Mock(
            constant_prefix=lambda: "test", constant_suffix=lambda: "-something.txt"
        ),
        mocker.Mock(
            constant_prefix=lambda: "test", constant_suffix=lambda: "-something.csv"
        ),
    ]
    return mocker.Mock(products=lambda: mock_products_rule_0, name="rule_0")


@pytest.fixture()
def rule_1(mocker):
    mock_products_rule_1 = [
        mocker.Mock(
            constant_prefix=lambda: "other/file-",
            constant_suffix=lambda: "-something.csv",
        )
    ]
    return mocker.Mock(products=lambda: mock_products_rule_1, name="rule_1")


@pytest.fixture()
def rule_2(mocker):
    mock_products_rule_2 = [
        mocker.Mock(
            constant_prefix=lambda: "other/file-dupe-",
            constant_suffix=lambda: "-something.csv",
        )
    ]
    return mocker.Mock(products=lambda: mock_products_rule_2, name="rule_2")


@pytest.fixture()
def rule_empty_suffix(mocker):
    return mocker.Mock(
        products=lambda: [
            Mock(constant_prefix=lambda: "test", constant_suffix=lambda: "")
        ]
    )


@pytest.fixture()
def rule_empty_prefix(mocker):
    return mocker.Mock(
        products=lambda: [
            Mock(constant_prefix=lambda: "", constant_suffix=lambda: "something.txt")
        ]
    )


@pytest.fixture()
def output_index(rule_0, rule_1, rule_2):
    return OutputIndex(rules=[rule_0, rule_1, rule_2])


@pytest.mark.parametrize(
    "target,expected_rules",
    [
        ("test-{wildcard}-something.txt", {"rule_0"}),
        ("test-{wildcard}-something.csv", {"rule_0"}),
        ("test-something.txt", {"rule_0"}),
        ("other/file--something.csv", {"rule_1"}),
        ("other/file-something.csv", {"rule_1"}),
        ("other/not-file-something.csv", set()),
        ("other/file-dupe-something.csv", {"rule_1", "rule_2"}),
    ],
)
def test_match(request, output_index, target, expected_rules):
    expected = {request.getfixturevalue(rule) for rule in expected_rules}
    assert output_index.match(target) == expected


@pytest.mark.parametrize(
    "target,expected_rules",
    [
        ("test-something", {"rule_empty_suffix"}),
        ("something", set()),
        ("something.txt", {"rule_empty_prefix"}),
        ("test-something.txt", {"rule_empty_suffix", "rule_empty_prefix"}),
    ],
)
def test_match_with_empty_components(
    request, rule_empty_suffix, rule_empty_prefix, target, expected_rules
):
    output_index = OutputIndex([rule_empty_suffix, rule_empty_prefix])
    expected = {request.getfixturevalue(rule) for rule in expected_rules}
    assert output_index.match(target) == expected


def test_empty_pattern_matches_everything(mocker):
    """Test that empty patterns match any filename"""
    rule = mocker.Mock(
        products=lambda: [Mock(constant_prefix=lambda: "", constant_suffix=lambda: "")]
    )
    output_index = OutputIndex([rule])
    assert rule in output_index.match("")
    assert rule in output_index.match("anything")
    assert rule in output_index.match("file.txt")


def test_empty_prefix_and_suffix(mocker):
    """Test with empty prefix and suffix"""
    rule = mocker.Mock(
        products=lambda: [Mock(constant_prefix=lambda: "", constant_suffix=lambda: "")]
    )
    output_index = OutputIndex([rule])
    matches = output_index.match("anything.txt")
    assert rule in matches


def test_case_sensitivity(mocker):
    """Test case sensitivity of matching"""
    rule = mocker.Mock(
        products=lambda: [
            Mock(constant_prefix=lambda: "Test", constant_suffix=lambda: "TXT")
        ]
    )
    output_index = OutputIndex([rule])

    matches_exact = output_index.match("Test.TXT")
    matches_lower = output_index.match("test.txt")
    assert rule in matches_exact
    assert rule not in matches_lower


@pytest.mark.parametrize(
    "target,expected_match",
    [
        ("test.txt", True),
        ("testing.txt", True),
        ("test.doc", False),
        ("other.txt", False),
        ("", False),
        ("test", False),
        ("bad.txt", False),
    ],
)
def test_parametrized_matches(mocker, target, expected_match):
    """Parametrized test for various matching scenarios"""
    rule = mocker.Mock(
        products=lambda: [
            Mock(constant_prefix=lambda: "test", constant_suffix=lambda: "txt")
        ]
    )
    output_index = OutputIndex([rule])
    matches = output_index.match(target)
    assert (rule in matches) == expected_match
