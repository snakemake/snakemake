import tokens
import textwrap

dd = textwrap.dedent

NEWLINE = set((tokens.NEWLINE, tokens.NL))

def newline(token, newline_tokens=set((tokens.NEWLINE, tokens.NL))):
    return token.type in newline_tokens


def indent(token):
    return token.type == tokens.indent


def dedent(token):
    return token.type == tokens.dedent


class StopAutomaton(Exception):
    pass

class TokenAutomaton:

    def __init__(self, tokenizer):
        self.tokenizer = tokenizer

    def consume(self):
        try:
            for token in self.tokenizer:
                for t in self.state(token):
                    yield t, token
        except StopAutomaton:
            return

    def error(self, msg, token):
        raise SyntaxError(msg,
            (self.tokenizer.path, token.start[0], None, None))


class RuleKeywordState(TokenAutomaton):

    def __init__(self, tokenizer):
        super().__init__(tokenizer)
        self.indent = 0
        self.line = 0

    @property
    def keyword(self):
        return self.__name__.lower()

    def start(self):
        yield "@workflow.{keyword}(".format(keyword=self.keyword)

    def end(self):
        yield ")"

    def colon(self, token):
        if token.string == ":":
            self.state = self.block
            for t in self.start():
                yield t
        else:
            self.error(
                "Colon expected after rule keyword {}.".format(self.keyword),
                token)

    def block(self, token):
        if newline(token):
            self.line += 1
            yield token.string
        elif indent(token):
            self.indent += 1
            yield token.string
        elif dedent(token):
            self.indent -= 1
            yield token.string
            if self.indent == 0:
                for t in self.end():
                    yield t
                raise StopAutomaton()
        else:
            yield token.string

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

    def block(self, token):
        for t in super().block(token):
            yield t
        self.shellcmd.append(token.string)
