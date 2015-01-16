__author__ = "Johannes Köster"


from collections import defaultdict


class Node:
    __slots__ = ["rules", "children"]
    def __init__(self):
        self.rules = set()
        self.children = defaultdict(Node)


class OutputIndex:
    def __init__(self, rules):
        self.root = Node()

        for rule in rules:
            output = list(rule.output)
            if rule.benchmark:
                output.append(rule.benchmark)
            for f in output:
                self.add_output(rule, f)

    def add_output(self, rule, f):
        node = self.root
        for c in f.constant_prefix():
            node = node.children[c]
            if rule in node.rules:
                # a prefix of file f is already recorded for this rule
                # hence we can stop here
                return
        node.rules.add(rule)

    def match(self, f):
        node = self.root
        for c in f:
            for rule in node.rules:
                yield rule
            node = node.children.get(c, None)
            if node is None:
                return
        for rule in node.rules:
            yield rule
