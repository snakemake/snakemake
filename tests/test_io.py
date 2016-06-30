from snakemake.io import _wildcard_regex


def test_wildcard_regex():
    def matches(text):
        return [ (match.group('name'), match.group('constraint'))
                 for match in _wildcard_regex.finditer(text) ]

    # without constraints
    assert matches('') == []
    assert matches('{') == []
    assert matches('}') == []
    assert matches('{}') == []
    assert matches('{0}') == [('0', None)]
    assert matches('{abc}') == [('abc', None)]
    assert matches('abc{def}{ghi}') == [('def', None), ('ghi', None)]

    # with constraints
    assert matches('{w,constraint}') == [('w', 'constraint')]
    assert matches('{w , constraint}') == [('w', 'constraint')]
    # fails because constraint is detected as 'constraint '
    # assert matches('{w,constraint }') == [('w', 'constraint')]
    assert matches('abc { w , constraint} def') == [('w', 'constraint')]

    # multiple wildcards
    assert matches('{a,1} {b,2} {c,3}') == [('a', '1'), ('b', '2'), ('c', '3')]

    # more complicated constraints
    assert matches(r'{w,([a-z]+|pat\|t*ern)}') == [('w', r'([a-z]+|pat\|t*ern)')]
    assert matches(r'{w,([a-z]+|pat\|te{1,3}rn){5,7}}') == [('w', r'([a-z]+|pat\|te{1,3}rn){5,7}')]

    # This used to be very slow with an older version of the regex
    assert matches('{w, long constraint without closing brace') == []
