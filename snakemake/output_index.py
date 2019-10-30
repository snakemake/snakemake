__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018-2019, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"

from itertools import chain

from snakemake.io import _IOFile


class OutputIndex:
    def __init__(self, rules):
        import datrie

        def prefixes(rule):
            return (str(o.constant_prefix()) for o in rule.products)

        def reverse_suffixes(rule):
            return (str(o.constant_suffix())[::-1] for o in rule.products)

        def calc_trie(subpatterns):
            t = datrie.Trie("".join(p for rule in rules for p in subpatterns(rule)))
            empty = list()
            for rule in rules:
                has_empty = False
                for p in subpatterns(rule):
                    if not p:
                        has_empty = True
                    if p not in t:
                        t[p] = [rule]
                    else:
                        t[p].append(rule)
                if has_empty:
                    empty.append(rule)
            return t, empty

        self.prefix_trie, self.empty_prefix = calc_trie(prefixes)
        self.suffix_trie, self.empty_suffix = calc_trie(reverse_suffixes)

    def match(self, targetfile):
        def match_pattern(pattern, trie, empty):
            return chain(chain.from_iterable(trie.iter_prefix_values(pattern)), empty)

        f = str(targetfile)
        hits = set(match_pattern(f, self.prefix_trie, self.empty_prefix))
        return hits.intersection(
            match_pattern(f[::-1], self.suffix_trie, self.empty_suffix)
        )
