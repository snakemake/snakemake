__author__ = "Johannes Köster"
__contributors__ = ["Per Unneberg"]
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import json
import re
import inspect
import textwrap
from itertools import chain
from collections import Mapping
import multiprocessing
import string
import shlex
import sys

from snakemake.io import regex, Namedlist, Wildcards
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError
import snakemake


def simplify_path(path):
    """Return a simplified version of the given path."""
    relpath = os.path.relpath(path)
    if relpath.startswith("../../"):
        return path
    else:
        return relpath


def linecount(filename):
    """Return the number of lines of given file.

    Args:
        filename (str): the path to the file
    """
    with open(filename) as f:
        return sum(1 for l in f)


def listfiles(pattern, restriction=None, omit_value=None):
    """Yield a tuple of existing filepaths for the given pattern.

    Wildcard values are yielded as the second tuple item.

    Args:
        pattern (str):       a filepattern. Wildcards are specified in snakemake syntax, e.g. "{id}.txt"
        restriction (dict):  restrict to wildcard values given in this dictionary
        omit_value (str):    wildcard value to omit

    Yields:
        tuple: The next file matching the pattern, and the corresponding wildcards object
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    if first_wildcard:
        dirname = os.path.dirname(pattern[:first_wildcard.start()])
        if not dirname:
            dirname = "."
    else:
        dirname = os.path.dirname(pattern)
    pattern = re.compile(regex(pattern))
    for dirpath, dirnames, filenames in os.walk(dirname):
        for f in chain(filenames, dirnames):
            if dirpath != ".":
                f = os.path.normpath(os.path.join(dirpath, f))
            match = re.match(pattern, f)
            if match:
                wildcards = Namedlist(fromdict=match.groupdict())
                if restriction is not None:
                    invalid = any(omit_value not in v and v != wildcards[k]
                                  for k, v in restriction.items())
                    if not invalid:
                        yield f, wildcards
                else:
                    yield f, wildcards


def makedirs(dirnames):
    """Recursively create the given directory or directories without
    reporting errors if they are present.
    """
    if isinstance(dirnames, str):
        dirnames = [dirnames]
    for dirname in dirnames:
        os.makedirs(dirname, exist_ok=True)


def report(text, path,
           stylesheet=os.path.join(os.path.dirname(__file__), "report.css"),
           defaultenc="utf8",
           template=None,
           metadata=None, **files):
    """Create an HTML report using python docutils.

    Attention: This function needs Python docutils to be installed for the
    python installation you use with Snakemake.

    All keywords not listed below are intepreted as paths to files that shall
    be embedded into the document. They keywords will be available as link
    targets in the text. E.g. append a file as keyword arg via F1=input[0]
    and put a download link in the text like this:

    .. code:: python

        report('''
        ==============
        Report for ...
        ==============

        Some text. A link to an embedded file: F1_.

        Further text.
        ''', outputpath, F1=input[0])

        Instead of specifying each file as a keyword arg, you can also expand
        the input of your rule if it is completely named, e.g.:

        report('''
        Some text...
        ''', outputpath, **input)

    Args:
        text (str):         The "restructured text" as it is expected by python docutils.
        path (str):         The path to the desired output file
        stylesheet (str):   An optional path to a css file that defines the style of the document. This defaults to <your snakemake install>/report.css. Use the default to get a hint how to create your own.
        defaultenc (str):   The encoding that is reported to the browser for embedded text files, defaults to utf8.
        template (str):     An optional path to a docutils HTML template.
        metadata (str):     E.g. an optional author name or email address.

    """
    try:
        import snakemake.report
    except ImportError:
        raise WorkflowError(
            "Python 3 package docutils needs to be installed to use the report function.")
    snakemake.report.report(text, path,
                            stylesheet=stylesheet,
                            defaultenc=defaultenc,
                            template=template,
                            metadata=metadata, **files)


def R(code):
    """Execute R code

    This function executes the R code given as a string.
    The function requires rpy2 to be installed.

    Args:
        code (str): R code to be executed
    """
    try:
        import rpy2.robjects as robjects
    except ImportError:
        raise ValueError(
            "Python 3 package rpy2 needs to be installed to use the R function.")
    robjects.r(format(textwrap.dedent(code), stepout=2))


class SequenceFormatter(string.Formatter):
    """string.Formatter subclass with special behavior for sequences.

    This class delegates formatting of individual elements to another
    formatter object. Non-list objects are formatted by calling the
    delegate formatter's "format_field" method. List-like objects
    (list, tuple, set, frozenset) are formatted by formatting each
    element of the list according to the specified format spec using
    the delegate formatter and then joining the resulting strings with
    a separator (space by default).

    """
    def __init__(self, separator=" ", element_formatter=string.Formatter(),
                 *args, **kwargs):
        self.separator = separator
        self.element_formatter = element_formatter

    def format_element(self, elem, format_spec):
        """Format a single element

        For sequences, this is called once for each element in a
        sequence. For anything else, it is called on the entire
        object. It is intended to be overridden in subclases.

        """
        return self.element_formatter.format_field(elem, format_spec)

    def format_field(self, value, format_spec):
        if isinstance(value, Wildcards):
            return ",".join("{}={}".format(name, value.get(name)) for name in value.keys())
        if isinstance(value, (list, tuple, set, frozenset)):
            return self.separator.join(self.format_element(v, format_spec) for v in value)
        else:
            return self.format_element(value, format_spec)


class QuotedFormatter(string.Formatter):
    """Subclass of string.Formatter that supports quoting.

    Using this formatter, any field can be quoted after formatting by
    appending "q" to its format string. By default, shell quoting is
    performed using "shlex.quote", but you can pass a different
    quote_func to the constructor. The quote_func simply has to take a
    string argument and return a new string representing the quoted
    form of the input string.

    Note that if an element after formatting is the empty string, it
    will not be quoted.

    """

    def __init__(self, quote_func=shlex.quote, *args, **kwargs):
        self.quote_func = quote_func
        super().__init__(*args, **kwargs)

    def format_field(self, value, format_spec):
        do_quote = format_spec.endswith("q")
        if do_quote:
            format_spec = format_spec[:-1]
        formatted = super().format_field(value, format_spec)
        if do_quote and formatted != '':
            formatted = self.quote_func(formatted)
        return formatted


class AlwaysQuotedFormatter(QuotedFormatter):
    """Subclass of QuotedFormatter that always quotes.

    Usage is identical to QuotedFormatter, except that it *always*
    acts like "q" was appended to the format spec.

    """

    def format_field(self, value, format_spec):
        if not format_spec.endswith("q"):
            format_spec += "q"
        return super().format_field(value, format_spec)


def format(_pattern, *args, stepout=1, _quote_all=False, **kwargs):
    """Format a pattern in Snakemake style.

    This means that keywords embedded in braces are replaced by any variable
    values that are available in the current namespace.
    """

    frame = inspect.currentframe().f_back
    while stepout > 1:
        if not frame.f_back:
            break
        frame = frame.f_back
        stepout -= 1

    variables = dict(frame.f_globals)
    # add local variables from calling rule/function
    variables.update(frame.f_locals)
    if "self" in variables and sys.version_info < (3, 5):
        # self is the first arg of fmt.format as well. Not removing it would
        # cause a multiple values error on Python <=3.4.2.
        del variables["self"]
    variables.update(kwargs)
    fmt = SequenceFormatter(separator=" ")
    if _quote_all:
        fmt.element_formatter = AlwaysQuotedFormatter()
    else:
        fmt.element_formatter = QuotedFormatter()
    try:
        return fmt.format(_pattern, *args, **variables)
    except KeyError as ex:
        raise NameError("The name {} is unknown in this context. Please "
                        "make sure that you defined that variable. "
                        "Also note that braces not used for variable access "
                        "have to be escaped by repeating them, "
                        "i.e. {{{{print $1}}}}".format(str(ex)))


class Unformattable:
    def __init__(self, errormsg="This cannot be used for formatting"):
        self.errormsg = errormsg

    def __str__(self):
        raise ValueError(self.errormsg)


def read_job_properties(jobscript,
                        prefix="# properties",
                        pattern=re.compile("# properties = (.*)")):
    """Read the job properties defined in a snakemake jobscript.

    This function is a helper for writing custom wrappers for the
    snakemake --cluster functionality. Applying this function to a
    jobscript will return a dict containing information about the job.
    """
    with open(jobscript) as jobscript:
        for m in map(pattern.match, jobscript):
            if m:
                return json.loads(m.group(1))


def min_version(version):
    """Require minimum snakemake version, raise workflow error if not met."""
    import pkg_resources
    if pkg_resources.parse_version(
        snakemake.__version__) < pkg_resources.parse_version(version):
        raise WorkflowError(
            "Expecting Snakemake version {} or higher.".format(version))


def update_config(config, overwrite_config):
    """Recursively update dictionary config with overwrite_config.

    See
    http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for details.

    Args:
      config (dict): dictionary to update
      overwrite_config (dict): dictionary whose items will overwrite those in config

    """

    def _update(d, u):
        for (key, value) in u.items():
            if (isinstance(value, Mapping)):
                d[key] = _update(d.get(key, {}), value)
            else:
                d[key] = value
        return d

    _update(config, overwrite_config)


def set_temporary_output(*rules):
    """Set the output of rules to temporary"""
    for rule in rules:
        logger.debug(
            "setting output of rule '{rule}' to temporary".format(rule=rule))
        rule.temp_output = set(rule.output)


def set_protected_output(*rules):
    """Set the output of rules to protected"""
    for rule in rules:
        logger.debug(
            "setting output of rule '{rule}' to protected".format(rule=rule))
        rule.protected_output = set(rule.output)


def available_cpu_count():
    """
    Return the number of available virtual or physical CPUs on this system.
    The number of available CPUs can be smaller than the total number of CPUs
    when the cpuset(7) mechanism is in use, as is the case on some cluster
    systems.

    Adapted from http://stackoverflow.com/a/1006301/715090
    """
    try:
        with open('/proc/self/status') as f:
            status = f.read()
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', status)
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return min(res, multiprocessing.cpu_count())
    except IOError:
        pass

    return multiprocessing.cpu_count()
