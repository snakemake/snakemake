__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import traceback
from tokenize import TokenError

from snakemake.logging import logger


def format_error(ex, lineno,
                 linemaps=None,
                 snakefile=None,
                 show_traceback=False):
    if linemaps is None:
        linemaps = dict()
    msg = str(ex)
    if linemaps and snakefile and snakefile in linemaps:
        lineno = linemaps[snakefile][lineno]
        if isinstance(ex, SyntaxError):
            msg = ex.msg
    location = (" in line {} of {}".format(lineno, snakefile) if
                lineno and snakefile else "")
    tb = ""
    if show_traceback:
        tb = "\n".join(format_traceback(cut_traceback(ex), linemaps=linemaps))
    return '{}{}{}{}'.format(ex.__class__.__name__, location, ":\n" + msg
                             if msg else ".", "\n{}".format(tb) if
                             show_traceback and tb else "")


def get_exception_origin(ex, linemaps):
    for file, lineno, _, _ in reversed(traceback.extract_tb(ex.__traceback__)):
        if file in linemaps:
            return lineno, file


def cut_traceback(ex):
    snakemake_path = os.path.dirname(__file__)
    for line in traceback.extract_tb(ex.__traceback__):
        dir = os.path.dirname(line[0])
        if not dir:
            dir = "."
        if not os.path.isdir(dir) or not os.path.samefile(snakemake_path, dir):
            yield line


def format_traceback(tb, linemaps):
    for file, lineno, function, code in tb:
        if file in linemaps:
            lineno = linemaps[file][lineno]
        if code is not None:
            yield '  File "{}", line {}, in {}'.format(file, lineno, function)


def log_verbose_traceback(ex):
    tb = "Full " + "".join(traceback.format_exception(type(ex), ex, ex.__traceback__))
    logger.debug(tb)


def print_exception(ex, linemaps):
    """
    Print an error message for a given exception.

    Arguments
    ex -- the exception
    linemaps -- a dict of a dict that maps for each snakefile
        the compiled lines to source code lines in the snakefile.
    """
    log_verbose_traceback(ex)
    if isinstance(ex, SyntaxError) or isinstance(ex, IndentationError):
        logger.error(format_error(ex, ex.lineno,
                                  linemaps=linemaps,
                                  snakefile=ex.filename,
                                  show_traceback=True))
        return
    origin = get_exception_origin(ex, linemaps)
    if origin is not None:
        lineno, file = origin
        logger.error(format_error(ex, lineno,
                                  linemaps=linemaps,
                                  snakefile=file,
                                  show_traceback=True))
        return
    elif isinstance(ex, TokenError):
        logger.error(format_error(ex, None, show_traceback=False))
    elif isinstance(ex, MissingRuleException):
        logger.error(format_error(ex, None,
                                  linemaps=linemaps,
                                  snakefile=ex.filename,
                                  show_traceback=False))
    elif isinstance(ex, RuleException):
        for e in ex._include + [ex]:
            if not e.omit:
                logger.error(format_error(e, e.lineno,
                                          linemaps=linemaps,
                                          snakefile=e.filename,
                                          show_traceback=True))
    elif isinstance(ex, WorkflowError):
        logger.error(format_error(ex, ex.lineno,
                                  linemaps=linemaps,
                                  snakefile=ex.snakefile,
                                  show_traceback=True))
    elif isinstance(ex, KeyboardInterrupt):
        logger.info("Cancelling snakemake on user request.")
    else:
        traceback.print_exception(type(ex), ex, ex.__traceback__)


class WorkflowError(Exception):
    @staticmethod
    def format_arg(arg):
        if isinstance(arg, str):
            return arg
        else:
            return "{}: {}".format(arg.__class__.__name__, str(arg))

    def __init__(self, *args, lineno=None, snakefile=None, rule=None):
        super().__init__("\n".join(self.format_arg(arg) for arg in args))
        if rule is not None:
            self.lineno = rule.lineno
            self.snakefile = rule.snakefile
        else:
            self.lineno = lineno
            self.snakefile = snakefile
        self.rule = rule


class WildcardError(WorkflowError):
    pass


class RuleException(Exception):
    """
    Base class for exception occuring withing the
    execution or definition of rules.
    """

    def __init__(self,
                 message=None,
                 include=None,
                 lineno=None,
                 snakefile=None,
                 rule=None):
        """
        Creates a new instance of RuleException.

        Arguments
        message -- the exception message
        include -- iterable of other exceptions to be included
        lineno -- the line the exception originates
        snakefile -- the file the exception originates
        """
        super(RuleException, self).__init__(message)
        self._include = set()
        if include:
            for ex in include:
                self._include.add(ex)
                self._include.update(ex._include)
        if rule is not None:
            if lineno is None:
                lineno = rule.lineno
            if snakefile is None:
                snakefile = rule.snakefile

        self._include = list(self._include)
        self.lineno = lineno
        self.filename = snakefile
        self.omit = not message

    @property
    def messages(self):
        return map(str, (ex for ex in self._include + [self] if not ex.omit))


class InputFunctionException(WorkflowError):
    def __init__(self, msg, wildcards=None, lineno=None, snakefile=None, rule=None):
        msg = self.format_arg(msg) + "\nWildcards:\n" + "\n".join(
            "{}={}".format(name, value) for name, value in wildcards.items())
        super().__init__(msg, lineno=lineno, snakefile=snakefile, rule=rule)


class MissingOutputException(RuleException):
    pass


class IOException(RuleException):
    def __init__(self, prefix, rule, files,
                 include=None,
                 lineno=None,
                 snakefile=None):
        message = ("{} for rule {}:\n{}".format(prefix, rule, "\n".join(files))
                   if files else "")
        super().__init__(message=message,
                         include=include,
                         lineno=lineno,
                         snakefile=snakefile,
                         rule=rule)


class MissingInputException(IOException):
    def __init__(self, rule, files, include=None, lineno=None, snakefile=None):
        super().__init__("Missing input files", rule, files, include,
                         lineno=lineno,
                         snakefile=snakefile)


class PeriodicWildcardError(RuleException):
    pass


class ProtectedOutputException(IOException):
    def __init__(self, rule, files, include=None, lineno=None, snakefile=None):
        super().__init__("Write-protected output files", rule, files, include,
                         lineno=lineno,
                         snakefile=snakefile)


class UnexpectedOutputException(IOException):
    def __init__(self, rule, files, include=None, lineno=None, snakefile=None):
        super().__init__("Unexpectedly present output files "
                         "(accidentally created by other rule?)", rule, files,
                         include,
                         lineno=lineno,
                         snakefile=snakefile)

class ImproperShadowException(RuleException):
    def __init__(self, rule, lineno=None, snakefile=None):
        super().__init__("Rule cannot shadow if using ThreadPoolExecutor",
                         rule=rule, lineno=lineno, snakefile=snakefile)


class AmbiguousRuleException(RuleException):
    def __init__(self, filename, job_a, job_b, lineno=None, snakefile=None):
        super().__init__(
            "Rules {job_a} and {job_b} are ambiguous for the file {f}.\n"
            "Expected input files:\n"
            "\t{job_a}: {job_a.input}\n"
            "\t{job_b}: {job_b.input}".format(job_a=job_a,
                                              job_b=job_b,
                                              f=filename),
            lineno=lineno,
            snakefile=snakefile)
        self.rule1, self.rule2 = job_a.rule, job_b.rule


class CyclicGraphException(RuleException):
    def __init__(self, repeatedrule, file, rule=None):
        super().__init__("Cyclic dependency on rule {}.".format(repeatedrule),
                         rule=rule)
        self.file = file


class MissingRuleException(RuleException):
    def __init__(self, file, lineno=None, snakefile=None):
        super().__init__(
            "No rule to produce {} (if you use input functions make sure that they don't raise unexpected exceptions).".format(
                file),
            lineno=lineno,
            snakefile=snakefile)


class UnknownRuleException(RuleException):
    def __init__(self, name, prefix="", lineno=None, snakefile=None):
        msg = "There is no rule named {}.".format(name)
        if prefix:
            msg = "{} {}".format(prefix, msg)
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class NoRulesException(RuleException):
    def __init__(self, lineno=None, snakefile=None):
        super().__init__("There has to be at least one rule.",
                         lineno=lineno,
                         snakefile=snakefile)


class IncompleteFilesException(RuleException):
    def __init__(self, files):
        super().__init__(
            "The files below seem to be incomplete. "
            "If you are sure that certain files are not incomplete, "
            "mark them as complete with\n\n"
            "    snakemake --cleanup-metadata <filenames>\n\n"
            "To re-generate the files rerun your command with the "
            "--rerun-incomplete flag.\nIncomplete files:\n{}".format(
                "\n".join(files)))


class IOFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)

class RemoteFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)

class HTTPFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)

class FTPFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)

class S3FileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)

class SFTPFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)

class DropboxFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)

class ClusterJobException(RuleException):
    def __init__(self, job, jobid, jobscript):
        super().__init__(
            "Error executing rule {} on cluster (jobid: {}, jobscript: {}). "
            "For detailed error see the cluster log.".format(job.rule.name,
                                                             jobid, jobscript),
            lineno=job.rule.lineno,
            snakefile=job.rule.snakefile)


class CreateRuleException(RuleException):
    pass


class TerminatedException(Exception):
    pass
