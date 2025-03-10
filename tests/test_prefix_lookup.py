from snakemake.common.prefix_lookup import PrefixLookup


def test_basic_match():
    lookup = PrefixLookup([("a", 1), ("b", 2)])
    assert lookup.match("a") == {1}
    assert lookup.match("b") == {2}


def test_multiple_matches():
    lookup = PrefixLookup([("", 0), ("a", 1), ("ab", 2), ("abc", 3)])
    assert lookup.match("abcd") == {0, 1, 2, 3}
    assert lookup.match("abc") == {0, 1, 2, 3}
    assert lookup.match("ab") == {0, 1, 2}
    assert lookup.match("a") == {0, 1}
    assert lookup.match("") == {0}


def test_no_matches():
    lookup = PrefixLookup([("a", 1), ("b", 2)])
    assert lookup.match("c") == set()
    assert lookup.match("cc") == set()


def test_empty_lookup():
    lookup = PrefixLookup([])
    assert lookup.match("abc") == set()


def test_overlapping_prefixes():
    lookup = PrefixLookup(
        [
            ("a", 1),
            ("aa", 2),
            ("a", 3),
        ]
    )
    assert lookup.match("aaa") == {1, 2, 3}
    assert lookup.match("a") == {1, 3}


def test_case_sensitivity():
    lookup = PrefixLookup([("a", 1), ("A", 2), ("ab", 3)])
    assert lookup.match("a") == {1}
    assert lookup.match("A") == {2}
    assert lookup.match("ab") == {3, 1}


def test_special_characters():
    lookup = PrefixLookup([(" ", 1), ("\t", 2), ("\n", 3), ("$#@!", 4), ("test$", 5)])
    assert lookup.match(" hello") == {1}
    assert lookup.match("\tworld") == {2}
    assert lookup.match("\ntest") == {3}
    assert lookup.match("$#@!test") == {4}
    assert lookup.match("test$abc") == {5}


def test_unicode_characters():
    lookup = PrefixLookup([("é", 1), ("ñ", 2), ("🌟", 3), ("é日本", 4)])
    assert lookup.match("étest") == {1}
    assert lookup.match("ñtest") == {2}
    assert lookup.match("🌟test") == {3}
    assert lookup.match("é日本語") == {1, 4}


def test_sorting_behavior():
    # Test that internal sorting doesn't affect matching
    lookup1 = PrefixLookup([("b", 1), ("a", 2), ("c", 3)])
    lookup2 = PrefixLookup([("a", 2), ("c", 3), ("b", 1)])
    assert lookup1.match("b") == lookup2.match("b")
    assert lookup1.match("abc") == lookup2.match("abc")


def test_different_value_types():
    lookup = PrefixLookup([("a", 42), ("b", "hello"), ("c", (1, 2, 3)), ("d", None)])
    assert lookup.match("abc") == {42}
    assert lookup.match("bcd") == {"hello"}
    assert lookup.match("cde") == {(1, 2, 3)}
    assert lookup.match("def") == {None}


def test_edge_cases():
    lookup = PrefixLookup(
        [
            ("", 1),
            (" ", 2),
            ("  ", 3),
            ("   ", 4),
        ]
    )
    assert lookup.match("abc") == {1}
    assert lookup.match(" abc") == {1, 2}
    assert lookup.match("  abc") == {1, 2, 3}
    assert lookup.match("   abc") == {1, 2, 3, 4}
    assert lookup.match("    abc") == {1, 2, 3, 4}


def test_prefix_matching_non_consecutive():
    """Test that demonstrates why we need to keep checking even after seeing a longer prefix,
    when prefixes aren't consecutive (like 'ab' vs 'abbc')."""
    lookup = PrefixLookup([("a", 1), ("ab", 2), ("abbcc", 3), ("abbc", 4), ("abbd", 5)])

    # Query "abbz" process:
    # 1. Try "abbd" (len 4) - no match, last_len = 4
    # 2. Try "abbc" (len 4) - no match, last_len = 4
    # 3. Try "ab" (len 2) - MATCHES, last_len = 2
    # 4. Try "a" (len 1) - MATCHES, last_len = 1
    result = lookup.match("abbz")
    assert result == {1, 2}  # Matches "a" and "ab"

    assert lookup.match("abbc") == {1, 2, 4}  # Matches "a", "ab", and "abbc"
    assert lookup.match("abbd") == {1, 2, 5}  # Matches "a", "ab", and "abbd"
    assert lookup.match("abbcc") == {
        1,
        2,
        3,
        4,
    }  # Matches "a", "ab", "abbcc", and "abbc"


def test_first_match_is_empty_string():
    lookup = PrefixLookup([("", 1)])
    matches = lookup.match("anything.txt")
    assert matches == {1}
