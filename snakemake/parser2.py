# -*- coding: utf-8 -*-

import tokenize
import textwrap
from itertools import chain


__author__ = "Johannes Köster"


dd = textwrap.dedent


def newline(token, newline_tokens=set((tokenize.NEWLINE, tokenize.NL))):
    return token.type in newline_tokens


def indent(token):
    return token.type == tokenize.INDENT


def dedent(token):
    return token.type == tokenize.DEDENT


def greater(token):
    return token.type == tokenize.OP and token.string == ">"


def name(token):
    return token.type == tokenize.NAME


def colon(token):
    return token.string == ":"


def comment(token):
    return token.type == tokenize.COMMENT


def string(token):
    return token.type == tokenize.STRING


def lineno(token):
    return token.start[0]


class StopAutomaton(Exception):

    def __init__(self, token):
        self.token = token


class TokenAutomaton:

    def __init__(self, tokenizer):
        self.tokenizer = tokenizer
        self.state = None

    def consume(self):
        for token in self.tokenizer:
            for t in self.state(token):
                yield t

    def error(self, msg, token):
        raise SyntaxError(msg,
            (self.tokenizer.path, lineno(token), None, None))


class KeywordState(TokenAutomaton):

    def __init__(self, tokenizer):
        super().__init__(tokenizer)
        self.indent = 0
        self.line = 0
        self.state = self.colon

    @property
    def keyword(self):
        return self.__class__.__name__.lower()

    def end(self):
        yield ")"

    def colon(self, token):
        if colon(token):
            self.state = self.block
            for t in self.start():
                yield t, token
        else:
            self.error(
                "Colon expected after keyword {}.".format(self.keyword),
                token)

    def block(self, token):
        if newline(token):
            self.line += 1
            yield token.string, token
        elif indent(token):
            self.indent += 1
            yield token.string, token
        elif dedent(token):
            self.indent -= 1
            yield token.string, token
        elif self.indent == 0 and self.line != 0:
                for t in self.end():
                    yield t, token
                raise StopAutomaton(token)
        else:
            for t in self.block_content(token):
                yield t

    def block_content(self, token):
        yield token.string, token


class GlobalKeywordState(KeywordState):

    def start(self):
        yield "workflow.{keyword}(".format(keyword=self.keyword)


class RuleKeywordState(KeywordState):

    def __init__(self, tokenizer, rulename=None):
        super().__init__(tokenizer)
        self.rulename = rulename

    def start(self):
        yield "@workflow.{keyword}(".format(keyword=self.keyword)


# Global keyword states


class Include(GlobalKeywordState):
    pass


class Workdir(GlobalKeywordState):
    pass


class Ruleorder(GlobalKeywordState):

    def block_content(self, token):
        if greater(token):
            yield ",", token
        elif name(token):
            yield token.string, token
        else:
            self.error('Expected a descending order of rule names, '
                'e.g. rule1 > rule2 > rule3 ...', token)


# Rule keyword states


class Input(RuleKeywordState):
    pass


class Output(RuleKeywordState):
    pass


class Params(RuleKeywordState):
    pass


class Threads(RuleKeywordState):
    pass


class Priority(RuleKeywordState):
    pass


class Version(RuleKeywordState):
    pass


class Log(RuleKeywordState):
    pass


class Message(RuleKeywordState):
    pass


class Run(RuleKeywordState):
    def __init__(self, tokenizer, rulename):
        super().__init__(tokenizer)
        self.rulename = rulename

    def start(self):
        yield textwrap.dedent("""
        @workflow.run
        __{rulename}(
        """).format(rulename=self.rulename)


class Shell(Run):

    def __init__(self, tokenizer, rulename):
        super().__init__(tokenizer, rulename)
        self.shellcmd = list()

    def start(self):
        yield "@workflow.shellcmd("

    def end(self):
        yield ")\n"
        for t in chain(
            super().start(),
            self.shellcmd,
            super().end()):
            yield t

    def block_content(self, token):
        self.shellcmd.append(token.string)
        yield token.string, token


class Rule(GlobalKeywordState):
    keywords = dict(
        input=Input,
        output=Output,
        params=Params,
        threads=Threads,
        priority=Priority,
        version=Version,
        log=Log,
        message=Message,
        run=Run,
        shell=Shell)

    def __init__(self, tokenizer):
        super().__init__(tokenizer)
        self.state = self.name

    def start(self):
        yield ("@workflow.rule(name={rulename}, lineno={lineno}, "
            "snakefile={snakefile})".format(
                rulename=self.rulename,
                lineno=self.lineno,
                snakefile=self.tokenizer.path))

    def end(self):
        yield ""

    def name(self, token):
        if name(token):
            self.rulename = token.string
        elif colon(token):
            self.lineno = lineno(token)
            self.state = self.block
            for t in self.start():
                yield t, token
        else:
            self.error("Expected name or colon after rule keyword.", token)

    def block_content(self, token):
        #import pdb; pdb.set_trace()
        if name(token):
            try:
                subautomaton = self.keywords[token.string](
                    self.tokenizer, rulename=self.rulename)
                for t in subautomaton.consume():
                    yield t
            except KeyError:
                self.error("Unexpected keyword {} in "
                    "rule definition".format(token.string), token)
            except StopAutomaton as e:
                # finished with keyword, go on in this state
                for t in self.block(e.token):
                    yield t
        elif comment(token):
            yield token.string, token
        elif string(token):
            yield "@workflow.docstring({})".format(token.string), token
        else:
            self.error("Expecting rule keyword, comment or docstrings "
                "inside a rule definition.", token)


class Python(TokenAutomaton):

    keywords = dict(
        include=Include,
        workdir=Workdir,
        ruleorder=Ruleorder,
        rule=Rule)

    def __init__(self, tokenizer):
        super().__init__(tokenizer)
        self.state = self.python

    def python(self, token):
        try:
            subautomaton = self.keywords[token.string](self.tokenizer)
            for t in subautomaton.consume():
                yield t
        except KeyError:
            yield token.string, token
        except StopAutomaton as e:
            for t in self.python(e.token):
                yield t


class Tokenizer:

    def __init__(self, path):
        self.path = path
        self.file = open(self.path)
        self.tokens = tokenize.generate_tokens(self.file.readline)

    def __next__(self):
        return next(self.tokens)

    def __iter__(self):
        return self

    def __exit__(self):
        self.file.close()


def parse(path):
    tokenizer = Tokenizer(path)
    automaton = Python(tokenizer)
    linemap = dict()
    compilation = list()
    lines = 1
    for t, orig_token in automaton.consume():
        l = lineno(orig_token)
        linemap.update(
            dict((i, l) for i in range(lines, lines + t.count("\n"))))
        lines += t.count("\n")
        compilation.append(t)
    print("".join(compilation))  # TODO implement correct spacing of tokens!
    return "".join(compilation), linemap
