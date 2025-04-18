__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import traceback
from tokenize import TokenError
from snakemake_interface_common.exceptions import WorkflowError, ApiError
from snakemake_interface_logger_plugins.common import LogEvent


def format_error(
    ex, lineno, linemaps=None, snakefile=None, show_traceback=False, rule=None
):
    if linemaps is None:
        linemaps = dict()
    msg = str(ex)
    if linemaps and snakefile and snakefile in linemaps:
        lineno = linemaps[snakefile][lineno]

    if isinstance(ex, SyntaxError):
        msg = ex.msg.split("(")[0]
        msg = f"{msg}:\n{ex.text}"

    location = ""
    if lineno and snakefile:
        location = f' in file "{snakefile}", line {lineno}'
        if rule:
            location = f" in rule {rule}{location}"

    tb = ""
    if show_traceback:
        tb = "\n".join(format_traceback(cut_traceback(ex), linemaps=linemaps))
    return "{}{}{}{}".format(
        ex.__class__.__name__,
        location,
        ":\n" + msg if msg else ".",
        f"\n{tb}" if show_traceback and tb else "",
    )


def get_exception_origin(ex, linemaps):
    for file, lineno, _, _ in reversed(traceback.extract_tb(ex.__traceback__)):
        if file in linemaps:
            return lineno, file


def cut_traceback(ex):
    lines = []
    snakemake_path = os.path.dirname(__file__)
    not_seen_snakemake = True
    for line in traceback.extract_tb(ex.__traceback__)[::-1]:
        dir = os.path.dirname(line[0])
        if not dir:
            dir = "."
        is_snakemake_dir = lambda path: os.path.realpath(path).startswith(
            os.path.realpath(snakemake_path)
        )
        if is_snakemake_dir(dir):
            not_seen_snakemake = False
        if not os.path.isdir(dir) or not_seen_snakemake:
            lines.append(line)
    return lines[::-1]


def format_traceback(tb, linemaps):
    for file, lineno, function, code in tb:
        if file in linemaps:
            lineno = linemaps[file][lineno]
        if code is not None:
            yield f'  File "{file}", line {lineno}, in {function}'


def log_verbose_traceback(ex):
    from snakemake.logging import logger

    tb = "Full " + "".join(traceback.format_exception(type(ex), ex, ex.__traceback__))
    logger.debug(tb)


def format_exception_to_string(ex, linemaps=None):
    """
    Returns the error message for a given exception as a string.

    Arguments
    ex -- the exception
    linemaps -- a dict of a dict that maps for each snakefile
        the compiled lines to source code lines in the snakefile.
    """
    if isinstance(ex, SyntaxError) or isinstance(ex, IndentationError):
        return format_error(
            ex,
            ex.lineno,
            linemaps=linemaps,
            snakefile=ex.filename,
            show_traceback=True,
        )

    origin = get_exception_origin(ex, linemaps) if linemaps is not None else None
    if origin is not None:
        lineno, file = origin
        return format_error(
            ex, lineno, linemaps=linemaps, snakefile=file, show_traceback=True
        )
    elif isinstance(ex, TokenError):
        return format_error(ex, None, show_traceback=False)
    elif isinstance(ex, MissingRuleException):
        return format_error(
            ex, None, linemaps=linemaps, snakefile=ex.filename, show_traceback=False
        )
    elif isinstance(ex, RuleException):
        error_string = ""
        for e in ex._include:
            if not e.omit:
                error_string += (
                    format_error(
                        e,
                        e.lineno,
                        linemaps=linemaps,
                        snakefile=e.filename,
                        show_traceback=True,
                    )
                    + "\n"
                )
        error_string += format_error(
            ex,
            ex.lineno,
            linemaps=linemaps,
            snakefile=ex.filename,
            show_traceback=True,
            rule=ex.rule,
        )
        return error_string
    elif isinstance(ex, WorkflowError):
        return format_error(
            ex,
            ex.lineno,
            linemaps=linemaps,
            snakefile=ex.snakefile,
            show_traceback=True,
            rule=ex.rule,
        )
    elif isinstance(ex, ApiError):
        return f"Error: {ex}"
    elif isinstance(ex, CliException):
        return f"Error: {ex}"
    elif isinstance(ex, KeyboardInterrupt):
        return "Cancelling snakemake on user request."
    else:
        return "\n".join(traceback.format_exception(ex))


def print_exception_warning(ex, linemaps=None, footer_message=""):
    """
    Print an error message for a given exception using logger warning.

    Arguments
    ex -- the exception
    linemaps -- a dict of a dict that maps for each snakefile
        the compiled lines to source code lines in the snakefile.
    """
    from snakemake.logging import logger

    log_verbose_traceback(ex)
    logger.warning(f"{format_exception_to_string(ex, linemaps)}\n{footer_message}")


def print_exception(ex, linemaps=None):
    """
    Print an error message for a given exception.

    Arguments
    ex -- the exception
    linemaps -- a dict of a dict that maps for each snakefile
        the compiled lines to source code lines in the snakefile.
    """
    from snakemake.logging import logger

    log_verbose_traceback(ex)
    logger.error(
        format_exception_to_string(ex, linemaps),
        extra=dict(event=LogEvent.ERROR, exception=ex.__class__.__name__),
    )


def update_lineno(ex: SyntaxError, linemaps):
    if ex.filename and ex.lineno:
        linemap = linemaps[ex.filename]
        try:
            ex.lineno = linemap[ex.lineno]
        except KeyError:
            # linemap does not yet contain the line, it must happen during parsing
            # such that no update is needed.
            pass
        return ex


class SourceFileError(WorkflowError):
    def __init__(self, msg):
        super().__init__(f"Error in source file definition: {msg}")


class WildcardError(WorkflowError):
    pass


class RuleException(Exception):
    """
    Base class for exception occurring within the
    execution or definition of rules.
    """

    def __init__(
        self, message=None, include=None, lineno=None, snakefile=None, rule=None
    ):
        """
        Creates a new instance of RuleException.

        Arguments
        message -- the exception message
        include -- iterable of other exceptions to be included
        lineno -- the line the exception originates
        snakefile -- the file the exception originates
        """
        super(RuleException, self).__init__(message)
        _include = set()
        if include:
            for ex in include:
                _include.add(ex)
                _include.update(ex._include)
        if rule is not None:
            if lineno is None:
                lineno = rule.lineno
            if snakefile is None:
                snakefile = rule.snakefile

        self._include = list(_include)
        self.rule = rule
        self.lineno = lineno
        self.filename = snakefile
        self.omit = not message

    @property
    def messages(self):
        return map(str, (ex for ex in self._include + [self] if not ex.omit))


class InputFunctionException(WorkflowError):
    def __init__(self, msg, wildcards=None, lineno=None, snakefile=None, rule=None):
        fmt_msg = (
            "Error:\n  "
            + self.format_arg(msg)
            + "\nWildcards:\n"
            + "\n".join(f"  {name}={value}" for name, value in wildcards.items())
        )
        if isinstance(msg, Exception):
            fmt_msg += "\nTraceback:\n" + "\n".join(
                format_traceback(cut_traceback(msg), rule.workflow.linemaps)
            )
        super().__init__(fmt_msg, lineno=lineno, snakefile=snakefile, rule=rule)


class ChildIOException(WorkflowError):
    def __init__(
        self,
        parent=None,
        child=None,
        wildcards=None,
        lineno=None,
        snakefile=None,
        rule=None,
    ):
        msg = "File/directory is a child to another output:\n" + "{}\n{}".format(
            parent, child
        )
        super().__init__(msg, lineno=lineno, snakefile=snakefile, rule=rule)


class IOException(RuleException):
    def __init__(self, prefix, job, files, include=None, lineno=None, snakefile=None):
        from snakemake.logging import format_wildcards
        from snakemake.io.fmt import fmt_iofile

        msg = ""
        if files:
            msg = f"{prefix} for rule {job.rule}:"
            if job.output:
                msg += "\n" + f"    output: {', '.join(map(fmt_iofile, job.output))}"
            if job.wildcards:
                msg += "\n" + f"    wildcards: {format_wildcards(job.wildcards)}"
            msg += "\n    affected files:\n        " + "\n        ".join(
                map(fmt_iofile, files)
            )
        super().__init__(
            message=msg,
            include=include,
            lineno=lineno,
            snakefile=snakefile,
            rule=job.rule,
        )


class MissingOutputException(RuleException):
    def __init__(
        self,
        message=None,
        include=None,
        lineno=None,
        snakefile=None,
        rule=None,
        jobid="",
    ):
        if jobid:
            jobid = f"{jobid} "
        message = f"Job {jobid} completed successfully, but some output files are missing. {message}"
        super().__init__(message, include, lineno, snakefile, rule)


class MissingInputException(IOException):
    def __init__(self, job, files, include=None, lineno=None, snakefile=None):
        msg = "Missing input files"

        if any(map(lambda f: f.startswith("~"), files)):
            msg += (
                "(Using '~' in your paths is not allowed as such platform "
                "specific syntax is not resolved by Snakemake. In general, "
                "try sticking to relative paths for everything inside the "
                "working directory.)"
            )
        super().__init__(msg, job, files, include, lineno=lineno, snakefile=snakefile)


class PeriodicWildcardError(RuleException):
    pass


class ProtectedOutputException(IOException):
    def __init__(self, job, files, include=None, lineno=None, snakefile=None):
        super().__init__(
            "Write-protected output files",
            job,
            files,
            include,
            lineno=lineno,
            snakefile=snakefile,
        )


class ImproperOutputException(IOException):
    def __init__(self, job, files, include=None, lineno=None, snakefile=None):
        super().__init__(
            "Outputs of incorrect type (directories when expecting files or vice versa). "
            "Output directories must be flagged with directory().",
            job,
            files,
            include,
            lineno=lineno,
            snakefile=snakefile,
        )


class ImproperShadowException(RuleException):
    def __init__(self, rule, lineno=None, snakefile=None):
        super().__init__(
            "Rule cannot shadow if using ThreadPoolExecutor",
            rule=rule,
            lineno=lineno,
            snakefile=snakefile,
        )


class AmbiguousRuleException(RuleException):
    def __init__(self, filename, job_a, job_b, lineno=None, snakefile=None):
        from snakemake import utils
        from snakemake.io.fmt import fmt_iofile

        filename = fmt_iofile(filename)

        wildcards_a = utils.format("{}", job_a._format_wildcards)
        wildcards_b = utils.format("{}", job_b._format_wildcards)
        super().__init__(
            "Rules {job_a} and {job_b} are ambiguous for the file {f}.\n"
            "Consider starting rule output with a unique prefix, constrain "
            "your wildcards, or use the ruleorder directive.\n"
            "Wildcards:\n"
            "\t{job_a}: {wildcards_a}\n"
            "\t{job_b}: {wildcards_b}\n"
            "Expected input files:\n"
            "\t{job_a}: {job_a.input}\n"
            "\t{job_b}: {job_b.input}\n"
            "Expected output files:\n"
            "\t{job_a}: {job_a.output}\n"
            "\t{job_b}: {job_b.output}".format(
                job_a=job_a,
                job_b=job_b,
                f=filename,
                wildcards_a=wildcards_a,
                wildcards_b=wildcards_b,
            ),
            lineno=lineno,
            snakefile=snakefile,
        )
        self.rule1, self.rule2 = job_a.rule, job_b.rule


class CyclicGraphException(RuleException):
    def __init__(self, repeatedrule, file, rule=None):
        super().__init__(f"Cyclic dependency on rule {repeatedrule}.", rule=rule)
        self.file = file


class MissingRuleException(RuleException):
    def __init__(self, file, lineno=None, snakefile=None):
        super().__init__(
            f"No rule to produce {file} (if you use input "
            "functions make sure that they don't raise unexpected exceptions).",
            lineno=lineno,
            snakefile=snakefile,
        )


class UnknownRuleException(RuleException):
    def __init__(self, name, prefix="", lineno=None, snakefile=None):
        msg = f"There is no rule named {name}."
        if prefix:
            msg = f"{prefix} {msg}"
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class NoRulesException(RuleException):
    def __init__(self, lineno=None, snakefile=None):
        super().__init__(
            "There has to be at least one rule.", lineno=lineno, snakefile=snakefile
        )


class IncompleteFilesException(RuleException):
    def __init__(self, files):
        from snakemake.io.fmt import fmt_iofile

        super().__init__(
            "The files below seem to be incomplete. "
            "If you are sure that certain files are not incomplete, "
            "mark them as complete with\n\n"
            "    snakemake --cleanup-metadata <filenames>\n\n"
            "To re-generate the files rerun your command with the "
            "--rerun-incomplete flag.\nIncomplete files:\n{}".format(
                "\n".join(map(fmt_iofile, files))
            )
        )


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


class AzureFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class SFTPFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class DropboxFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class XRootDFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class NCBIFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class WebDAVFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class ZenodoFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


class ClusterJobException(RuleException):
    def __init__(self, job_info, jobid):
        super().__init__(
            "Error executing rule {} on cluster (jobid: {}, external: {}, jobscript: {}). "
            "For detailed error see the cluster log.".format(
                job_info.job.rule.name, jobid, job_info.jobid, job_info.jobscript
            ),
            lineno=job_info.job.rule.lineno,
            snakefile=job_info.job.rule.snakefile,
        )


class CreateRuleException(RuleException):
    pass


class TerminatedException(Exception):
    pass


class CreateCondaEnvironmentException(WorkflowError):
    pass


class SpawnedJobError(Exception):
    pass


class CheckSumMismatchException(WorkflowError):
    """ "should be called to indicate that checksum of a file compared to known
    hash does not match, typically done with large downloads, etc.
    """

    pass


class IncompleteCheckpointException(Exception):
    def __init__(self, rule, targetfile):
        super().__init__(
            "The requested checkpoint output is not yet created. "
            "If you see this error, you have likely tried to use "
            "checkpoint output outside of an input function, or "
            "you have tried to call an input function directly "
            "via <function_name>(). Please check the docs at "
            "https://snakemake.readthedocs.io/en/stable/"
            "snakefiles/rules.html#data-dependent-conditional-execution "
            "and note that the input function in the example rule "
            "'aggregate' is NOT called, but passed to the rule "
            "by name, such that Snakemake can call it internally "
            "once the checkpoint is finished."
        )
        self.rule = rule
        from snakemake.io import checkpoint_target

        self.targetfile = checkpoint_target(targetfile)


class InputOpenException(Exception):
    def __init__(self, iofile):
        self.iofile = iofile
        self.rule = None


class CacheMissException(Exception):
    pass


class LockException(WorkflowError):
    def __init__(self):
        super().__init__(
            "Error: Directory cannot be locked. Please make "
            "sure that no other Snakemake process is trying to create "
            "the same files in the following directory:\n{}\n"
            "If you are sure that no other "
            "instances of snakemake are running on this directory, "
            "the remaining lock was likely caused by a kill signal or "
            "a power loss. It can be removed with "
            "the --unlock argument.".format(os.getcwd())
        )


class ResourceScopesException(Exception):
    def __init__(self, msg, invalid_resources):
        super().__init__(msg, invalid_resources)
        self.msg = msg
        self.invalid_resources = invalid_resources


class CliException(Exception):
    def __init__(self, msg):
        super().__init__(msg)
        self.msg = msg


class LookupError(WorkflowError):
    def __init__(self, msg=None, exc=None, query=None, dpath=None):
        msg = f" {msg}" if msg is not None else ""
        expr = ""
        if query is not None:
            expr = f" with query: {repr(query)}"
        if dpath is not None:
            expr = f" with dpath: {repr(dpath)}"
        annotated_msg = f"Error in lookup function{expr}.{msg}"
        args = [annotated_msg]
        if exc is not None:
            args.append(exc)
        super().__init__(*args)


class MissingOutputFileCachePathException(Exception):
    pass


def is_file_not_found_error(exc, considered_files):
    # TODO find a better way to detect whether the input files are not present
    if isinstance(exc, FileNotFoundError) and exc.filename in considered_files:
        return True
    elif isinstance(exc, WorkflowError) and "FileNotFoundError" in str(exc):
        return True
    else:
        return False
