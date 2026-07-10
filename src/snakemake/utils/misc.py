import os
import json
import re
import textwrap
from itertools import chain
import collections
import multiprocessing
import sys

from snakemake.iocontainers import Namedlist
from snakemake.common import ON_WINDOWS
from snakemake.exceptions import WorkflowError
from snakemake.utils.format import format


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
    from snakemake.io import regex_from_filepattern

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
    except ModuleNotFoundError:
        raise ValueError(
            "Python 3 package rpy2 needs to be installed to use the R function."
        )
    robjects.r(format(textwrap.dedent(code), stepout=2))


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
