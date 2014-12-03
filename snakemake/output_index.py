

from collections import defaultdict


class Node:
    __slots__ = ["rules", "children"]
    def __init__(self):
        self.rules = []
        self.children = defaultdict(Node)


class OutputIndex:
    def __init__(self, rules):
        self.root = Node()

        for rule in rules:
            for f in rule.output:
                node = self.root
                for c in f.constant_prefix():
                    node = node.children[c]
                node.rules.append(rule)

    def match(self, file):
        node = self.root
        for c in file.constant_prefix():
            for rule in node.rules:
                yield rule
            node = node.children.get(c, None)
            if node is None:
                return
