from itertools import chain

import datrie

from snakemake.io import _IOFile


class OutputIndex:
    def __init__(self, rules):
        def prefixes(rule):
            return map(_IOFile.constant_prefix, rule.products)
        def reverse_suffixes(rule):
            return (_IOFile.constant_suffix(o)[::-1] for o in rule.products)
        def calc_trie(subpatterns):
            t = datrie.Trie("".join(p for rule in rules
                                    for p in subpatterns(rule)))
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
            return chain(chain.from_iterable(trie.iter_prefix_values(pattern)),
                         empty)

        hits = set(match_pattern(str(targetfile),
                   self.prefix_trie, self.empty_prefix))
        hits = hits.intersection(match_pattern(targetfile[::-1],
                                               self.suffix_trie,
                                               self.empty_suffix))
        return hits
