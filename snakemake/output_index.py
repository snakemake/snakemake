__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

from collections import defaultdict

from snakemake.io import _IOFile


class Node:
    __slots__ = ["rules", "children"]

    def __init__(self):
        self.rules = set()
        self.children = defaultdict(Node)

    def __repr__(self):
        return "({}) -> {}".format(list(map(str, self.rules)), dict(self.children))


class OutputIndex:
    def __init__(self, rules):
        self.root = Node()

        for rule in rules:
            output = list(rule.output)
            if rule.benchmark:
                output.append(rule.benchmark)
            for constant_prefix in sorted(map(_IOFile.constant_prefix,
                                              output)):
                self.add_output(rule, constant_prefix)

    def add_output(self, rule, constant_prefix):
        node = self.root
        for c in constant_prefix:
            node = node.children[c]
            if rule in node.rules:
                # a prefix of file f is already recorded for this rule
                # hence we can stop here
                return
        node.rules.add(rule)

    def match(self, f):
        node = self.root
        rules = set()
        for c in f:
            for rule in node.rules:
                rules.add(rule)
            node = node.children.get(c, None)
            if node is None:
                return rules
        for rule in node.rules:
            rules.add(rule)
        return rules
