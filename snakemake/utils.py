# -*- coding: utf-8 -*-

__author__ = "Johannes Köster"

import os
import io
import re
import fnmatch
import mimetypes
import base64
import inspect
import textwrap
import tempfile
import subprocess
import shutil
import mimetypes
import datetime
from itertools import chain

from snakemake.io import regex, Namedlist
from snakemake.logging import logger


def linecount(filename):
    """
    Return the number of lines of given file

    Arguments
    filename -- the path to the file
    """
    with open(filename) as f:
        return sum(1 for l in f)


def listfiles(pattern, restriction=None, omit_value=None):
    """
    Yield a tuple of existing filepaths for the given pattern.
    Wildcard values are yielded as the second tuple item.

    Arguments
    pattern -- a filepattern.
        Wildcards are specified in snakemake syntax, e.g. "{id}.txt"
    """
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
                f = os.path.join(dirpath, f)
            match = re.match(pattern, f)
            if match and len(match.group()) == len(f):
                wildcards = Namedlist(fromdict=match.groupdict())
                if restriction is not None:
                    invalid = any(
                        omit_value not in v and v != wildcards[k]
                        for k, v in restriction.items())
                    if not invalid:
                        yield f, wildcards
                else:
                    yield f, wildcards


def makedirs(dirnames):
    """
    Recursively create the given directory or directories without
    reporting errors if they are present.
    """
    if isinstance(dirnames, str):
        dirnames = [dirnames]
    for dirname in dirnames:
        if not os.path.exists(dirname):
            os.makedirs(dirname)


def report(
    text, path,
    stylesheet=os.path.join(os.path.dirname(__file__), "report.css"),
    defaultenc="utf8", template=None, **files):
    """
    Create an HTML report using python docutils.
    Attention: This function needs Python docutils to be installed for the
    python installation you use with Snakemake.

    Arguments
    text -- The "restructured text" as it is expected by python docutils.
    path -- The path to the desired output file
    stylesheet -- An optional path to a css file that defines the style of the
        document. This defaults to <your snakemake install>/report.css.
        Use the default to get a hint how to create your own.
    defaultenc -- The encoding that is reported to the browser for embedded
        text files, defaults to utf8.
    template -- An optional path to a docutils HTML template.

    All other keyword args are intepreted as paths to files that shall be
    embedded into the document. They keywords will be available as link
    targets in the text. E.g. append a file as keyword arg via F1=input[0]
    and put a download link in the text like this:

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

    """
    outmime, _ = mimetypes.guess_type(path)
    if outmime != "text/html":
        raise ValueError("Path to report output has to be an HTML file.")
    from docutils.core import publish_file
    definitions = textwrap.dedent("""
    .. role:: raw-html(raw)
       :format: html

    """)

    metadata = textwrap.dedent("""

    .. container::
       :name: metadata

       {metadata}

    """).format(metadata=datetime.date.today().isoformat())

    text = format(textwrap.dedent(text), stepout=2)

    attachments = [textwrap.dedent("""
        .. container::
           :name: attachments
           
        """)]
    for name, file in sorted(files.items()):
        mime, encoding = mimetypes.guess_type(file)
        if mime is None:
            mime = "text/plain"
            logger.warning("Could not detect mimetype for {}, assuming "
            "text/plain.".format(file))
        if encoding is None:
            encoding = defaultenc
        with open(file, "rb") as f:
            data = base64.b64encode(f.read())
        attachments.append(
            '''
   .. container::
      :name: {name}

      [{name}] :raw-html:`<a href="data:{mime};charset={charset};filename={filename};base64,{data}" download="{filename}" draggable="true">{filename}</a>`
            '''.format(
                name=name,
                filename=os.path.basename(file),
                mime=mime,
                charset=encoding,
                data=data.decode()))

    text = definitions + text + "\n\n" + "\n\n".join(attachments) + metadata

    overrides = dict()
    if template is not None:
        overrides["template"] = template
    if stylesheet is not None:
        overrides["stylesheet_path"] = stylesheet
    html = open(path, "w")
    publish_file(
        source=io.StringIO(text), destination=html,
        writer_name="html", settings_overrides=overrides)


def R(code):
    import rpy2.robjects as robjects
    robjects.r(format(textwrap.dedent(code), stepout=2))


def format(string, *args, stepout=1, **kwargs):
    class SequenceFormatter:
        def __init__(self, sequence):
            self._sequence = sequence

        def __getitem__(self, i):
            return self._sequence[i]

        def __str__(self):
            return " ".join(self._sequence)

    frame = inspect.currentframe().f_back
    while stepout > 1:
        if not frame.f_back:
            break
        frame = frame.f_back
        stepout -= 1

    variables = dict(frame.f_globals)
    # add local variables from calling rule/function
    variables.update(frame.f_locals)
    variables.update(kwargs)
    strmethods = list()
    for key, value in list(variables.items()):
        if type(value) in (list, tuple, set, frozenset):
            variables[key] = SequenceFormatter(value)
    try:
        return string.format(*args, **variables)
    except KeyError as ex:
        raise NameError(
            "The name {} is unknown in this context. Please"
            "make sure that you defined that variable. "
            "Also note that braces not used for variable access "
            "have to be escaped by repeating them, "
            "i.e. {{print $1}}".format(str(ex)))


class Unformattable:

    def __init__(self, errormsg="This cannot be used for formatting"):
        self.errormsg = errormsg

    def __str__(self):
        raise ValueError(self.errormsg)
