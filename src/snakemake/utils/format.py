# compat: python 3.7 (script support)
# Must be kept compatible to Python 3.7 because it is used in Snakemake's
# Python script support.
# Only modules from python standard library or other Snakemake modules that are
# Python 3.7 compatible should be imported here (except for methods
# that are only called by Snakemake itself)!

from snakemake.iocontainers import Wildcards
import string
import shlex
import inspect
import sys


def format(_pattern, *args, stepout=1, _quote_all=False, quote_func=None, **kwargs):
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
        fmt.element_formatter = AlwaysQuotedFormatter(quote_func)
    else:
        fmt.element_formatter = QuotedFormatter(quote_func)
    try:
        return fmt.format(_pattern, *args, **variables)
    except KeyError as ex:
        if (
            "wildcards" in variables
            and str(ex).strip("'") in variables["wildcards"].keys()
        ):
            raise NameError(
                "The name '{0}' is unknown in this context. "
                "Did you mean 'wildcards.{0}'?".format(str(ex).strip("'"))
            )
        raise NameError(
            "The name {} is unknown in this context. Please "
            "make sure that you defined that variable. "
            "Also note that braces not used for variable access "
            "have to be escaped by repeating them, "
            "i.e. {{{{print $1}}}}".format(str(ex))
        )


class Unformattable:
    def __init__(self, errormsg="This cannot be used for formatting"):
        self.errormsg = errormsg

    def __str__(self):
        raise ValueError(self.errormsg)


class SequenceFormatter(string.Formatter):
    """string.Formatter subclass with special behavior for sequences.

    This class delegates the formatting of individual elements to another
    formatter object. Non-list objects are formatted by calling the
    delegate formatter's "format_field" method. List-like objects
    (list, tuple, set, frozenset) are formatted by formatting each
    element of the list according to the specified format spec using
    the delegate formatter and then joining the resulting strings with
    a separator (space by default).

    """

    def __init__(
        self, separator=" ", element_formatter=string.Formatter(), *args, **kwargs
    ):
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
            return ",".join(
                f"{name}={value}"
                for name, value in sorted(value.items(), key=lambda item: item[0])
            )
        if isinstance(value, (list, tuple, set, frozenset)):
            return self.separator.join(
                self.format_element(v, format_spec) for v in value
            )
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

    def __init__(self, quote_func=None, *args, **kwargs):
        if quote_func is None:
            quote_func = shlex.quote
        self.quote_func = quote_func
        super().__init__(*args, **kwargs)

    def format_field(self, value, format_spec):
        if format_spec.endswith("u"):
            format_spec = format_spec[:-1]
            do_quote = False
        else:
            do_quote = format_spec.endswith("q")
            if do_quote:
                format_spec = format_spec[:-1]
        formatted = super().format_field(value, format_spec)
        if do_quote and formatted != "":
            formatted = self.quote_func(formatted)
        return formatted


class AlwaysQuotedFormatter(QuotedFormatter):
    """Subclass of QuotedFormatter that always quotes.

    Usage is identical to QuotedFormatter, except that it *always*
    acts like "q" was appended to the format spec, unless u (for unquoted) is appended.

    """

    def format_field(self, value, format_spec):
        if not format_spec.endswith("q") and not format_spec.endswith("u"):
            format_spec += "q"
        return super().format_field(value, format_spec)
