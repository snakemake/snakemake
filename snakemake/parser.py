__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import tokenize
import textwrap
import os
from urllib.error import HTTPError, URLError, ContentTooShortError
import urllib.request
from io import TextIOWrapper

from snakemake.exceptions import WorkflowError

dd = textwrap.dedent

INDENT = "\t"


def is_newline(token, newline_tokens=set((tokenize.NEWLINE, tokenize.NL))):
    return token.type in newline_tokens


def is_indent(token):
    return token.type == tokenize.INDENT


def is_dedent(token):
    return token.type == tokenize.DEDENT


def is_op(token):
    return token.type == tokenize.OP


def is_greater(token):
    return is_op(token) and token.string == ">"


def is_comma(token):
    return is_op(token) and token.string == ","


def is_name(token):
    return token.type == tokenize.NAME


def is_colon(token):
    return is_op(token) and token.string == ":"


def is_comment(token):
    return token.type == tokenize.COMMENT


def is_string(token):
    return token.type == tokenize.STRING


def is_eof(token):
    return token.type == tokenize.ENDMARKER


def lineno(token):
    return token.start[0]


class StopAutomaton(Exception):
    def __init__(self, token):
        self.token = token


class TokenAutomaton:

    subautomata = dict()

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        self.root = root
        self.snakefile = snakefile
        self.state = None
        self.base_indent = base_indent
        self.line = 0
        self.indent = 0
        self.was_indented = False
        self.lasttoken = None
        self._dedent = dedent

    @property
    def dedent(self):
        return self._dedent

    @property
    def effective_indent(self):
        return self.base_indent + self.indent - self.dedent

    def indentation(self, token):
        if is_indent(token) or is_dedent(token):
            self.indent = token.end[1] - self.base_indent
            self.was_indented |= self.indent > 0

    def consume(self):
        for token in self.snakefile:
            self.indentation(token)
            try:
                for t, orig in self.state(token):
                    if self.lasttoken == "\n" and not t.isspace():
                        yield INDENT * self.effective_indent, orig
                    yield t, orig
                    self.lasttoken = t
            except tokenize.TokenError as e:
                self.error(
                    str(e).split(",")[0].strip("()''"), token
                )  # TODO the inferred line number seems to be wrong sometimes

    def error(self, msg, token):
        raise SyntaxError(msg,
                          (self.snakefile.path, lineno(token), None, None))

    def subautomaton(self, automaton, *args, **kwargs):
        return self.subautomata[automaton](
            self.snakefile, *args,
            base_indent=self.base_indent + self.indent,
            dedent=self.dedent,
            root=False, **kwargs)


class KeywordState(TokenAutomaton):

    prefix = ""

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile,
                         base_indent=base_indent,
                         dedent=dedent,
                         root=root)
        self.line = 0
        self.state = self.colon

    @property
    def keyword(self):
        return self.__class__.__name__.lower()[len(self.prefix):]

    def end(self):
        yield ")"

    def decorate_end(self, token):
        for t in self.end():
            yield t, token

    def colon(self, token):
        if is_colon(token):
            self.state = self.block
            for t in self.start():
                yield t, token
        else:
            self.error("Colon expected after keyword {}.".format(self.keyword),
                       token)

    def is_block_end(self, token):
        return (self.line and self.indent <= 0) or is_eof(token)

    def block(self, token):
        if self.lasttoken == "\n" and is_comment(token):
            # ignore lines containing only comments
            self.line -= 1
        if self.is_block_end(token):
            for t, token_ in self.decorate_end(token):
                yield t, token_
            yield "\n", token
            raise StopAutomaton(token)

        if is_newline(token):
            self.line += 1
            yield token.string, token
        elif not (is_indent(token) or is_dedent(token)):
            if is_comment(token):
                yield token.string, token
            else:
                for t in self.block_content(token):
                    yield t

    def yield_indent(self, token):
        return token.string, token

    def block_content(self, token):
        yield token.string, token


class GlobalKeywordState(KeywordState):
    def start(self):
        yield "workflow.{keyword}(".format(keyword=self.keyword)


class DecoratorKeywordState(KeywordState):
    decorator = None
    args = list()

    def start(self):
        yield "@workflow.{}".format(self.decorator)
        yield "\n"
        yield "def __{}({}):".format(self.decorator, ", ".join(self.args))

    def end(self):
        yield ""


class RuleKeywordState(KeywordState):
    def __init__(self, snakefile,
                 base_indent=0,
                 dedent=0,
                 root=True,
                 rulename=None):
        super().__init__(snakefile,
                         base_indent=base_indent,
                         dedent=dedent,
                         root=root)
        self.rulename = rulename

    def start(self):
        yield "\n"
        yield "@workflow.{keyword}(".format(keyword=self.keyword)


class SubworkflowKeywordState(KeywordState):
    prefix = "Subworkflow"

    def start(self):
        yield ", {keyword}=".format(keyword=self.keyword)

    def end(self):
        # no end needed
        return list()

# Global keyword states


class Include(GlobalKeywordState):
    pass


class Workdir(GlobalKeywordState):
    pass


class Configfile(GlobalKeywordState):
    pass


class Ruleorder(GlobalKeywordState):
    def block_content(self, token):
        if is_greater(token):
            yield ",", token
        elif is_name(token):
            yield '"{}"'.format(token.string), token
        else:
            self.error('Expected a descending order of rule names, '
                       'e.g. rule1 > rule2 > rule3 ...', token)

# subworkflows


class SubworkflowSnakefile(SubworkflowKeywordState):
    pass


class SubworkflowWorkdir(SubworkflowKeywordState):
    pass


class Subworkflow(GlobalKeywordState):

    subautomata = dict(snakefile=SubworkflowSnakefile,
                       workdir=SubworkflowWorkdir)

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile,
                         base_indent=base_indent,
                         dedent=dedent,
                         root=root)
        self.state = self.name
        self.has_snakefile = False
        self.has_workdir = False
        self.has_name = False
        self.primary_token = None

    def end(self):
        if not (self.has_snakefile or self.has_workdir):
            self.error(
                "A subworkflow needs either a path to a Snakefile or to a workdir.",
                self.primary_token)
        yield ")"

    def name(self, token):
        if is_name(token):
            yield "workflow.subworkflow('{name}'".format(
                name=token.string), token
            self.has_name = True
        elif is_colon(token) and self.has_name:
            self.primary_token = token
            self.state = self.block
        else:
            self.error("Expected name after subworkflow keyword.", token)

    def block_content(self, token):
        if is_name(token):
            try:
                if token.string == "snakefile":
                    self.has_snakefile = True
                if token.string == "workdir":
                    self.has_workdir = True
                for t in self.subautomaton(token.string).consume():
                    yield t
            except KeyError:
                self.error("Unexpected keyword {} in "
                           "subworkflow definition".format(token.string),
                           token)
            except StopAutomaton as e:
                self.indentation(e.token)
                for t in self.block(e.token):
                    yield t
        elif is_comment(token):
            yield "\n", token
            yield token.string, token
        elif is_string(token):
            # ignore docstring
            pass
        else:
            self.error("Expecting subworkflow keyword, comment or docstrings "
                       "inside a subworkflow definition.", token)


class Localrules(GlobalKeywordState):
    def block_content(self, token):
        if is_comma(token):
            yield ",", token
        elif is_name(token):
            yield '"{}"'.format(token.string), token
        else:
            self.error('Expected a comma separated list of rules that shall '
                       'not be executed by the cluster command.', token)

# Rule keyword states


class Input(RuleKeywordState):
    pass


class Output(RuleKeywordState):
    pass


class Params(RuleKeywordState):
    pass


class Threads(RuleKeywordState):
    pass


class Shadow(RuleKeywordState):
    pass


class Resources(RuleKeywordState):
    pass


class Priority(RuleKeywordState):
    pass


class Version(RuleKeywordState):
    pass


class Log(RuleKeywordState):
    pass


class Message(RuleKeywordState):
    pass


class Benchmark(RuleKeywordState):
    pass


class Run(RuleKeywordState):
    def __init__(self, snakefile, rulename,
                 base_indent=0,
                 dedent=0,
                 root=True):
        super().__init__(snakefile,
                         base_indent=base_indent,
                         dedent=dedent,
                         root=root)
        self.rulename = rulename

    def start(self):
        yield "@workflow.run"
        yield "\n"
        yield ("def __rule_{rulename}(input, output, params, wildcards, threads, "
               "resources, log, version):".format(rulename=self.rulename if self.rulename is not None else self.snakefile.rulecount))

    def end(self):
        yield ""

    def is_block_end(self, token):
        return (self.line and self.was_indented and self.indent <= 0
                ) or is_eof(token)


class Shell(Run):

    overwrite_shellcmd = None

    def __init__(self, snakefile, rulename,
                 base_indent=0,
                 dedent=0,
                 root=True):
        super().__init__(snakefile, rulename,
                         base_indent=base_indent,
                         dedent=dedent,
                         root=root)
        self.shellcmd = list()
        self.token = None
        if self.overwrite_shellcmd is not None:
            self.block_content = self.overwrite_block_content

    def is_block_end(self, token):
        return (self.line and self.indent <= 0) or is_eof(token)

    def start(self):
        yield "@workflow.shellcmd("

    def end(self):
        # the end is detected. So we can savely reset the indent to zero here
        self.indent = 0
        yield "\n"
        yield ")"
        yield "\n"
        for t in super().start():
            yield t
        yield "\n"
        yield INDENT * (self.effective_indent + 1)
        yield "shell("
        yield "\n".join(self.shellcmd)
        yield "\n"
        yield ")"
        for t in super().end():
            yield t

    def decorate_end(self, token):
        if self.token is None:
            # no block after shell keyword
            self.error(
                "Shell command must be given as string after the shell keyword.",
                token)
        for t in self.end():
            yield t, self.token

    def block_content(self, token):
        self.token = token
        self.shellcmd.append(token.string)
        yield token.string, token

    def overwrite_block_content(self, token):
        if self.token is None:
            self.token = token
            shellcmd = '"{}"'.format(self.overwrite_shellcmd)
            self.shellcmd.append(shellcmd)
            yield shellcmd, token


class Script(Run):
    def __init__(self, snakefile, rulename,
                 base_indent=0,
                 dedent=0,
                 root=True):
        super().__init__(snakefile, rulename,
                         base_indent=base_indent,
                         dedent=dedent,
                         root=root)
        self.path = list()
        self.token = None

    def is_block_end(self, token):
        return (self.line and self.indent <= 0) or is_eof(token)

    def start(self):
        for t in super().start():
            yield t
        yield "\n"
        yield INDENT * (self.effective_indent + 1)
        yield "script("
        yield '"{}"'.format(
            os.path.abspath(os.path.dirname(self.snakefile.path)))
        yield ", "

    def end(self):
        # the end is detected. So we can savely reset the indent to zero here
        self.indent = 0
        yield ", input, output, params, wildcards, threads, resources, log, config"
        yield ")"
        for t in super().end():
            yield t

    def decorate_end(self, token):
        if self.token is None:
            # no block after script keyword
            self.error(
                "Script path must be given as string after the script keyword.",
                token)
        for t in self.end():
            yield t, self.token

    def block_content(self, token):
        self.token = token
        self.path.append(token.string)
        yield token.string, token


class Wrapper(Script):
    def start(self):
        for t in super(Script, self).start():
            yield t
        yield "\n"
        yield INDENT * (self.effective_indent + 1)
        yield 'wrapper('


class Rule(GlobalKeywordState):
    subautomata = dict(input=Input,
                       output=Output,
                       params=Params,
                       threads=Threads,
                       resources=Resources,
                       priority=Priority,
                       version=Version,
                       log=Log,
                       message=Message,
                       benchmark=Benchmark,
                       shadow=Shadow,
                       run=Run,
                       shell=Shell,
                       script=Script,
                       wrapper=Wrapper)

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile,
                         base_indent=base_indent,
                         dedent=dedent,
                         root=root)
        self.state = self.name
        self.lineno = None
        self.rulename = None
        self.run = False
        self.snakefile.rulecount += 1

    def start(self):
        yield ("@workflow.rule(name={rulename}, lineno={lineno}, "
               "snakefile='{snakefile}')".format(
                   rulename=("'{}'".format(self.rulename) if self.rulename is
                             not None else None),
                   lineno=self.lineno,
                   snakefile=self.snakefile.path.replace('\\', '\\\\')))

    def end(self):
        if not self.run:
            yield "@workflow.norun()"
            yield "\n"
            for t in self.subautomaton("run", rulename=self.rulename).start():
                yield t
            # the end is detected.
            # So we can savely reset the indent to zero here
            self.indent = 0
            yield "\n"
            yield INDENT * (self.effective_indent + 1)
            yield "pass"

    def name(self, token):
        if is_name(token):
            self.rulename = token.string
        elif is_colon(token):
            self.lineno = self.snakefile.lines + 1
            self.state = self.block
            for t in self.start():
                yield t, token
        else:
            self.error("Expected name or colon after rule keyword.", token)

    def block_content(self, token):
        if is_name(token):
            try:
                if token.string == "run" or token.string == "shell" or token.string == "script" or token.string == "wrapper":
                    if self.run:
                        raise self.error(
                            "Multiple run or shell keywords in rule {}.".format(
                                self.rulename), token)
                    self.run = True
                for t in self.subautomaton(token.string,
                                           rulename=self.rulename).consume():
                    yield t
            except KeyError:
                self.error("Unexpected keyword {} in "
                           "rule definition".format(token.string), token)
            except StopAutomaton as e:
                self.indentation(e.token)
                for t in self.block(e.token):
                    yield t
        elif is_comment(token):
            yield "\n", token
            yield token.string, token
        elif is_string(token):
            yield "\n", token
            yield "@workflow.docstring({})".format(token.string), token
        else:
            self.error("Expecting rule keyword, comment or docstrings "
                       "inside a rule definition.", token)

    @property
    def dedent(self):
        return self.indent


class OnSuccess(DecoratorKeywordState):
    decorator = "onsuccess"
    args = ["log"]


class OnError(DecoratorKeywordState):
    decorator = "onerror"
    args = ["log"]

class OnStart(DecoratorKeywordState):
    decorator = "onstart"
    args = ["log"]


class Python(TokenAutomaton):

    subautomata = dict(include=Include,
                       workdir=Workdir,
                       configfile=Configfile,
                       ruleorder=Ruleorder,
                       rule=Rule,
                       subworkflow=Subworkflow,
                       localrules=Localrules,
                       onsuccess=OnSuccess,
                       onerror=OnError,
                       onstart=OnStart)

    def __init__(self, snakefile, base_indent=0, dedent=0, root=True):
        super().__init__(snakefile,
                         base_indent=base_indent,
                         dedent=dedent,
                         root=root)
        self.state = self.python

    def python(self, token):
        if not (is_indent(token) or is_dedent(token)):
            if self.lasttoken is None or self.lasttoken.isspace():
                try:
                    for t in self.subautomaton(token.string).consume():
                        yield t
                except KeyError:
                    yield token.string, token
                except StopAutomaton as e:
                    self.indentation(e.token)
                    for t in self.python(e.token):
                        yield t
            else:
                yield token.string, token


class Snakefile:
    def __init__(self, path, rulecount=0):
        self.path = path
        try:
            self.file = open(self.path, encoding="utf-8")
        except FileNotFoundError as e:
            try:
                self.file = TextIOWrapper(urllib.request.urlopen(self.path),
                                          encoding="utf-8")
            except (HTTPError, URLError, ContentTooShortError, ValueError):
                raise WorkflowError("Failed to open {}.".format(path))

        self.tokens = tokenize.generate_tokens(self.file.readline)
        self.rulecount = rulecount
        self.lines = 0

    def __next__(self):
        return next(self.tokens)

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


def format_tokens(tokens):
    t_ = None
    for t in tokens:
        if t_ and not t.isspace() and not t_.isspace():
            yield " "
        yield t
        t_ = t


def parse(path, overwrite_shellcmd=None, rulecount=0):
    Shell.overwrite_shellcmd = overwrite_shellcmd
    with Snakefile(path, rulecount=rulecount) as snakefile:
        automaton = Python(snakefile)
        linemap = dict()
        compilation = list()
        for t, orig_token in automaton.consume():
            l = lineno(orig_token)
            linemap.update(dict((i, l) for i in range(
                snakefile.lines + 1, snakefile.lines + t.count("\n") + 1)))
            snakefile.lines += t.count("\n")
            compilation.append(t)
        compilation = "".join(format_tokens(compilation))
        if linemap:
            last = max(linemap)
            linemap[last + 1] = linemap[last]
        return compilation, linemap, snakefile.rulecount
