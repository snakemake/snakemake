__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import json
import re
import inspect
import textwrap
from itertools import chain
import collections
import multiprocessing
import string
import shlex
import sys

from snakemake.io import Namedlist, Wildcards
from snakemake.common.configfile import _load_configfile
from snakemake.logging import logger
from snakemake.common import ON_WINDOWS
from snakemake.exceptions import WorkflowError
from snakemake.io import regex_from_filepattern


def validate(data, schema, set_default=True):
    """Validate data with JSON schema at given path.

    Args:
        data (object): data to validate. Can be a config dict or a pandas data frame.
        schema (str): Path to JSON schema used for validation. The schema can also be
            in YAML format. If validating a pandas data frame, the schema has to
            describe a row record (i.e., a dict with column names as keys pointing
            to row values). See https://json-schema.org. The path is interpreted
            relative to the Snakefile when this function is called.
        set_default (bool): set default values defined in schema. See
            https://python-jsonschema.readthedocs.io/en/latest/faq/ for more
            information
    """
    frame = inspect.currentframe().f_back
    workflow = frame.f_globals.get("workflow")

    if workflow and (workflow.modifier.skip_validation or workflow.remote_exec):
        # skip if a corresponding modifier has been defined or if this is a
        # remote job. In the latter case, the schema has been already  validated by the
        # main process.
        return

    from snakemake.sourcecache import LocalSourceFile, infer_source_file
    import jsonschema
    from jsonschema import Draft202012Validator, validators
    from referencing import Registry, Resource

    schemafile = infer_source_file(schema)

    if isinstance(schemafile, LocalSourceFile) and not schemafile.isabs() and workflow:
        # if workflow object is not available this has not been started from a workflow
        schemafile = workflow.current_basedir.join(schemafile.get_path_or_uri())

    source = (
        workflow.sourcecache.open(schemafile)
        if workflow
        else schemafile.get_path_or_uri()
    )

    schema = _load_configfile(source, filetype="Schema")

    def retrieve_uri(uri):
        return Resource.from_contents(contents=_load_configfile(uri, filetype="Schema"))

    resource = Resource.from_contents(contents=schema)
    registry = Registry(retrieve=retrieve_uri).with_resource(
        uri=schemafile.get_path_or_uri(), resource=resource
    )
    Validator = Draft202012Validator(schema, registry=registry)

    # Taken from https://python-jsonschema.readthedocs.io/en/latest/faq/
    def extend_with_default(validator_class):
        validate_properties = validator_class.VALIDATORS["properties"]

        def set_defaults(validator, properties, instance, schema):
            for property, subschema in properties.items():
                if "default" in subschema:
                    instance.setdefault(property, subschema["default"])

            for error in validate_properties(
                validator,
                properties,
                instance,
                schema,
            ):
                yield error

        return validators.extend(
            validator_class,
            {"properties": set_defaults},
        )

    if Validator.META_SCHEMA["$schema"] != schema["$schema"]:
        logger.warning(
            f"No validator found for JSON Schema version identifier '{schema['$schema']}'"
        )
        logger.warning(
            f"Defaulting to validator for JSON Schema version '{Validator.META_SCHEMA['$schema']}'"
        )
        logger.warning("Note that schema file may not be validated correctly.")
    Defaultvalidator = extend_with_default(Validator)

    def _validate_record(record):
        if set_default:
            Defaultvalidator(schema, registry=registry).validate(record)
            return record
        else:
            Validator.validate(record)

    def _validate_pandas(data):
        try:
            import pandas as pd

            if isinstance(data, pd.DataFrame):
                logger.debug("Validating pandas DataFrame")

                recordlist = []
                for i, record in enumerate(data.to_dict("records")):
                    # Exclude NULL values
                    record = {k: v for k, v in record.items() if pd.notnull(v)}
                    try:
                        recordlist.append(_validate_record(record))
                    except jsonschema.exceptions.ValidationError as e:
                        raise WorkflowError(
                            f"Error validating row {i} of data frame.", e
                        )

                if set_default:
                    newdata = pd.DataFrame(recordlist, data.index)
                    # Add missing columns
                    newcol = newdata.columns[~newdata.columns.isin(data.columns)]
                    data[newcol] = None
                    # Fill in None values with values from newdata
                    data.update(newdata)

            else:
                return False
        except ImportError:
            return False
        return True

    def _validate_polars(data):
        try:
            import polars as pl

            if isinstance(data, pl.DataFrame):
                logger.debug("Validating polars DataFrame")

                recordlist = []
                for i, record in enumerate(data.iter_rows(named=True)):
                    # Exclude NULL values
                    record = {
                        k: v
                        for k, v in record.items()
                        if pl.Series(k, [v]).is_not_null().all()
                    }
                    try:
                        recordlist.append(_validate_record(record))
                    except jsonschema.exceptions.ValidationError as e:
                        raise WorkflowError(
                            f"Error validating row {i} of data frame.", e
                        )

                if set_default:
                    newdata = pl.DataFrame(recordlist)
                    # Add missing columns
                    newcol = [col for col in newdata.columns if col not in data.columns]
                    [
                        data.insert_column(
                            len(data.columns),
                            pl.lit(None, newdata[col].dtype).alias(col),
                        )
                        for col in newcol
                    ]
                    # Fill in None values with values from newdata
                    for i in range(data.shape[0]):
                        for j in range(data.shape[1]):
                            if data[i, j] is None:
                                data[i, j] = newdata[i, j]

            elif isinstance(data, pl.LazyFrame):
                # If a LazyFrame is being used, probably it is a large dataframe (so check only first 1000 records)
                logger.debug("Validating first 1000 rows of polars LazyFrame")

                recordlist = []
                for i, record in enumerate(
                    data.head(1000).collect().iter_rows(named=True)
                ):
                    # Exclude NULL values
                    record = {
                        k: v
                        for k, v in record.items()
                        if pl.Series(k, [v]).is_not_null().all()
                    }
                    try:
                        recordlist.append(_validate_record(record))
                    except jsonschema.exceptions.ValidationError as e:
                        raise WorkflowError(
                            f"Error validating row {i} of data frame.", e
                        )

                if set_default:
                    logger.warning("LazyFrame does not support setting default values.")

            else:
                return False
        except ImportError:
            return False
        return True

    if isinstance(data, dict):
        logger.debug("Validating dictionary")
        try:
            _validate_record(data)
        except jsonschema.exceptions.ValidationError as e:
            raise WorkflowError("Error validating config file.", e)
        logger.debug("Dictionary validated!")
    else:
        if _validate_pandas(data):
            logger.debug("Pandas dataframe validated!")
        elif _validate_polars(data):
            logger.debug("Polars dataframe validated!")


def simplify_path(path):
    """Return a simplified version of the given path."""
    relpath = os.path.relpath(path)
    if relpath.startswith("../../"):
        return path
    else:
        return relpath


def linecount(filename):
    """Return the number of lines of the given file.

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
        dirname = os.path.dirname(pattern[: first_wildcard.start()])
        if not dirname:
            dirname = "."
    else:
        dirname = os.path.dirname(pattern)
    pattern = re.compile(regex_from_filepattern(pattern))

    for dirpath, dirnames, filenames in os.walk(dirname):
        for f in chain(filenames, dirnames):
            if dirpath != ".":
                f = os.path.normpath(os.path.join(dirpath, f))
            match = re.match(pattern, f)
            if match:
                wildcards = Namedlist(fromdict=match.groupdict())
                if restriction is not None:
                    invalid = any(
                        omit_value not in v and v != wildcards[k]
                        for k, v in restriction.items()
                    )
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


def report(
    text,
    path,
    stylesheet=None,
    defaultenc="utf8",
    template=None,
    metadata=None,
    **files,
):
    """Create an HTML report using python docutils.

    This is deprecated in favor of the --report flag.

    Attention: This function needs Python docutils to be installed for the
    python installation you use with Snakemake.

    All keywords not listed below are interpreted as paths to files that shall
    be embedded into the document. The keywords will be available as link
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
        stylesheet (str):   An optional path to a CSS file that defines the style of the document. This defaults to <your snakemake install>/report.css. Use the default to get a hint on how to create your own.
        defaultenc (str):   The encoding that is reported to the browser for embedded text files, defaults to utf8.
        template (str):     An optional path to a docutils HTML template.
        metadata (str):     E.g. an optional author name or email address.

    """
    if stylesheet is None:
        os.path.join(os.path.dirname(__file__), "report.css")
    try:
        import snakemake.report
    except ImportError:
        raise WorkflowError(
            "Python 3 package docutils needs to be installed to use the report function."
        )
    snakemake.report.report(
        text,
        path,
        stylesheet=stylesheet,
        defaultenc=defaultenc,
        template=template,
        metadata=metadata,
        **files,
    )


def R(code):
    """Execute R code.

    This is deprecated in favor of the ``script`` directive.
    This function executes the R code given as a string.
    The function requires rpy2 to be installed.

    Args:
        code (str): R code to be executed
    """
    try:
        import rpy2.robjects as robjects
    except ImportError:
        raise ValueError(
            "Python 3 package rpy2 needs to be installed to use the R function."
        )
    robjects.r(format(textwrap.dedent(code), stepout=2))


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


def read_job_properties(
    jobscript, prefix="# properties", pattern=re.compile("# properties = (.*)")
):
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
    from packaging.version import parse
    from snakemake.common import __version__

    if parse(__version__) < parse(version):
        raise WorkflowError(
            "Expecting Snakemake version {} or higher (you are currently using {}).".format(
                version, __version__
            )
        )


def update_config(config, overwrite_config):
    """Recursively update dictionary config with overwrite_config in-place.

    See
    https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for details.

    Args:
      config (dict): dictionary to update
      overwrite_config (dict): dictionary whose items will overwrite those in config
    """

    def _update_config(config, overwrite_config):
        """Necessary as recursive calls require a return value,
        but `update_config()` has no return value.
        """
        for key, value in overwrite_config.items():
            if not isinstance(config, collections.abc.Mapping):
                # the config cannot be updated as it is no dict
                # -> just overwrite it with the new value
                config = {}
            if isinstance(value, collections.abc.Mapping):
                sub_config = config.get(key, {})
                config[key] = _update_config(sub_config, value)
            else:
                config[key] = value
        return config

    _update_config(config, overwrite_config)


def available_cpu_count():
    """
    Return the number of available virtual or physical CPUs on this system.
    The number of available CPUs can be smaller than the total number of CPUs
    when the cpuset(7) mechanism is in use, as is the case on some cluster
    systems.

    Adapted from https://stackoverflow.com/a/1006301/715090
    """
    try:
        with open("/proc/self/status") as f:
            status = f.read()
        m = re.search(r"(?m)^Cpus_allowed:\s*(.*)$", status)
        if m:
            res = bin(int(m.group(1).replace(",", ""), 16)).count("1")
            if res > 0:
                return min(res, multiprocessing.cpu_count())
    except IOError:
        pass

    return multiprocessing.cpu_count()


def argvquote(arg, force=True):
    """Returns an argument quoted in such a way that CommandLineToArgvW
    on Windows will return the argument string unchanged.
    This is the same thing Popen does when supplied with a list of arguments.
    Arguments in a command line should be separated by spaces; this
    function does not add these spaces. This implementation follows the
    suggestions outlined here:
    https://blogs.msdn.microsoft.com/twistylittlepassagesallalike/2011/04/23/everyone-quotes-command-line-arguments-the-wrong-way/
    """
    if not force and len(arg) != 0 and not any([c in arg for c in ' \t\n\v"']):
        return arg
    else:
        n_backslashes = 0
        cmdline = '"'
        for c in arg:
            if c == "\\":
                # first count the number of current backslashes
                n_backslashes += 1
                continue
            if c == '"':
                # Escape all backslashes and the following double quotation mark
                cmdline += (n_backslashes * 2 + 1) * "\\"
            else:
                # backslashes are not special here
                cmdline += n_backslashes * "\\"
            n_backslashes = 0
            cmdline += c
        # Escape all backslashes, but let the terminating
        # double quotation mark we add below be interpreted
        # as a metacharacter
        cmdline += +n_backslashes * 2 * "\\" + '"'
        return cmdline


def cmd_exe_quote(arg):
    """Quotes an argument in a cmd.exe compliant way."""
    arg = argvquote(arg)
    cmd_exe_metachars = '^()%!"<>&|'
    for char in cmd_exe_metachars:
        arg.replace(char, "^" + char)
    return arg


def os_sync():
    """Ensure flush to disk"""
    if not ON_WINDOWS:
        os.sync()


def find_bash_on_windows():
    """
    Find the path to a usable bash on windows.
    The first attempt is to look for a bash installed with a git conda package.
    Alternatively, try bash installed with 'Git for Windows'.
    """
    if not ON_WINDOWS:
        return None
    # First look for bash in git's conda package
    bashcmd = os.path.join(os.path.dirname(sys.executable), r"Library\bin\bash.exe")
    if not os.path.exists(bashcmd):
        # Otherwise try bash installed with "Git for Windows".
        import winreg

        try:
            key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, "SOFTWARE\\GitForWindows")
            gfwp, _ = winreg.QueryValueEx(key, "InstallPath")
            bashcmd = os.path.join(gfwp, "bin\\bash.exe")
        except FileNotFoundError:
            bashcmd = ""
    return bashcmd if os.path.exists(bashcmd) else None


class Paramspace:
    """A wrapper for pandas dataframes that provides helpers for using them as a parameter
    space in Snakemake.

    This is heavily inspired by @soumitrakp work on JUDI (https://github.com/ncbi/JUDI).

    By default, a directory structure with on folder level per parameter is created
    (e.g. column1~{column1}/column2~{column2}/***).

    The exact behavior can be tweaked with four parameters:

      - ``filename_params`` takes a list of column names of the passed dataframe.
        These names are used to build the filename (separated by '_') in the order
        in which they are passed.
        All remaining parameters will be used to generate a directory structure.
        Example for a data frame with four columns named column1 to column4:

        | ``Paramspace(df, filename_params=["column3", "column2"])`` ->
        | column1~{value1}/column4~{value4}/column3~{value3}_column2~{value2}

        If ``filename_params="*"``, all columns of the dataframe are encoded into
        the filename instead of parent directories.

      - ``param_sep`` takes a string which is used to join the column name and
        column value in the generated paths (Default: '~'). Example:

        | ``Paramspace(df, param_sep=":")`` ->
        | column1:{value1}/column2:{value2}/column3:{value3}/column4:{value4}

      - ``filename_sep`` takes a string which is used to join the parameter
        entries listed in ``filename_params`` in the generated paths
        (Default: '_'). Example:

        | ``Paramspace(df, filename_params="*", filename_sep="-")`` ->
        | column1~{value1}-column2~{value2}-column3~{value3}-column4~{value4}

      - ``single_wildcard`` takes a string which is used to replace the
        default behavior of using a wildcard for each column in the dataframe
        with a single wildcard that is used to encode all column values.
        The given string is the name of that wildcard. The value of the wildcard
        for individual instances of the paramspace is still controlled by above
        other arguments. The single_wildcard mechanism can be handy if you want
        to define a rule that shall be used for multiple paramspaces with different
        columns.
    """

    def __init__(
        self,
        dataframe,
        filename_params=None,
        param_sep="~",
        filename_sep="_",
        single_wildcard=None,
    ):
        self.dataframe = dataframe
        self.param_sep = param_sep
        self.filename_sep = filename_sep
        self.single_wildcard = single_wildcard
        if filename_params is None or not filename_params:
            # create a pattern of the form {}/{}/{} with one entry for each
            # column in the dataframe
            self.pattern = "/".join([r"{}"] * len(self.dataframe.columns))
            self.ordered_columns = self.dataframe.columns
        else:
            if isinstance(filename_params, str) and filename_params == "*":
                filename_params = dataframe.columns

            if any((param not in dataframe.columns for param in filename_params)):
                raise KeyError(
                    "One or more entries of filename_params are not valid column names for the param file."
                )
            elif len(set(filename_params)) != len(filename_params):
                raise ValueError("filename_params must be unique")
            # create a pattern of the form {}/{}_{} with one entry for each
            # column in the dataframe. The number of underscore-separated
            # fields is equal to the number filename_params
            self.pattern = "/".join(
                [r"{}"] * (len(self.dataframe.columns) - len(filename_params) + 1)
            )
            self.pattern = self.filename_sep.join(
                [self.pattern] + [r"{}"] * (len(filename_params) - 1)
            )
            self.ordered_columns = [
                param
                for param in self.dataframe.columns
                if param not in filename_params
            ]
            self.ordered_columns.extend(list(filename_params))
        self.dataframe = self.dataframe[self.ordered_columns]

    @property
    def wildcard_pattern(self):
        """Wildcard pattern over all columns of the underlying dataframe of the form
        column1~{column1}/column2~{column2}/*** or of the provided custom pattern.
        """
        if self.single_wildcard:
            return f"{{{self.single_wildcard}}}"
        else:
            return self.pattern.format(
                *map(
                    self.param_sep.join(("{0}", "{{{0}}}")).format, self.ordered_columns
                )
            )

    @property
    def instance_patterns(self):
        """Iterator over all instances of the parameter space (dataframe rows),
        formatted as file patterns of the form column1~{value1}/column2~{value2}/...
        or of the provided custom pattern.
        """
        import pandas as pd

        fmt_value = lambda value: "NA" if pd.isna(value) else value
        return (
            self.pattern.format(
                *(
                    self.param_sep.join(("{}", "{}")).format(name, fmt_value(value))
                    for name, value in row._asdict().items()
                )
            )
            for row in self.dataframe.itertuples(index=False)
        )

    def instance(self, wildcards):
        """Obtain instance (dataframe row) with the given wildcard values."""
        import pandas as pd
        from snakemake.io import regex_from_filepattern

        def convert_value_dtype(name, value):
            if self.dataframe.dtypes[name] == bool and value == "False":
                # handle problematic case when boolean False is returned as
                # boolean True because the string "False" is misinterpreted
                return False
            if value == "NA":
                return pd.NA
            else:
                return pd.Series([value]).astype(self.dataframe.dtypes[name])[0]

        if self.single_wildcard:
            wildcard_value = wildcards.get(self.single_wildcard)
            if wildcard_value is None:
                raise WorkflowError(
                    f"Error processing paramspace: wildcard {self.single_wildcard} is not used in rule."
                )

            pattern = self.pattern.format(
                *(
                    f"{name}{self.param_sep}{{{name}}}"
                    for name in self.dataframe.columns
                )
            )
            rexp = re.compile(regex_from_filepattern(pattern))
            match = rexp.match(wildcard_value)
            if not match:
                raise WorkflowError(
                    f"Error processing paramspace: wildcard {self.single_wildcard}={wildcards.get(self.single_wildcard)} does not match pattern {pattern}."
                )
            return {
                name: convert_value_dtype(name, value)
                for name, value in match.groupdict().items()
            }
        else:
            return {
                name: convert_value_dtype(name, value)
                for name, value in wildcards.items()
                if name in self.ordered_columns
            }

    def __getattr__(self, name):
        import pandas as pd

        ret = getattr(self.dataframe, name)
        if isinstance(ret, pd.DataFrame):
            return Paramspace(ret)
        return ret

    def __getitem__(self, key):
        import pandas as pd

        ret = self.dataframe[key]
        if isinstance(ret, pd.DataFrame):
            return Paramspace(ret)
        return ret
