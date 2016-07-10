__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import mimetypes
import base64
import textwrap
import datetime
import io

from docutils.parsers.rst.directives.images import Image, Figure
from docutils.parsers.rst import directives
from docutils.core import publish_file

from snakemake.utils import format
from snakemake.logging import logger


class EmbeddedMixin(object):
    """
    Replaces the URI of a directive with a base64-encoded version.

    Useful for embedding images/figures in reports.
    """

    def run(self):
        """
        Image.run() handles most of the
        """
        result = Image.run(self)
        reference = directives.uri(self.arguments[0])
        self.options['uri'] = data_uri(reference)
        return result

# Create (and register) new image:: and figure:: directives that use a base64
# data URI instead of pointing to a filename.


class EmbeddedImage(Image, EmbeddedMixin):
    pass


directives.register_directive('embeddedimage', EmbeddedImage)


class EmbeddedFigure(Figure, EmbeddedMixin):
    pass


directives.register_directive('embeddedfigure', EmbeddedFigure)


def data_uri(file, defaultenc="utf8"):
    """Craft a base64 data URI from file with proper encoding and mimetype."""
    mime, encoding = mimetypes.guess_type(file)
    if mime is None:
        mime = "text/plain"
        logger.info("Could not detect mimetype for {}, assuming "
                    "text/plain.".format(file))
    if encoding is None:
        encoding = defaultenc
    with open(file, "rb") as f:
        data = base64.b64encode(f.read())
    uri = ("data:{mime};charset={charset};filename={filename};base64,{data}"
           "".format(filename=os.path.basename(file),
                     mime=mime,
                     charset=encoding,
                     data=data.decode()))
    return uri


def report(text, path,
           stylesheet=os.path.join(os.path.dirname(__file__), "report.css"),
           defaultenc="utf8",
           template=None,
           metadata=None, **files):
    outmime, _ = mimetypes.guess_type(path)
    if outmime != "text/html":
        raise ValueError("Path to report output has to be an HTML file.")
    definitions = textwrap.dedent("""
    .. role:: raw-html(raw)
       :format: html

    """)

    metadata = textwrap.dedent("""

    .. container::
       :name: metadata

       {metadata}{date}

    """).format(metadata=metadata + " | " if metadata else "",
                date=datetime.date.today().isoformat())

    text = format(textwrap.dedent(text), stepout=3)

    attachments = [textwrap.dedent("""
        .. container::
           :name: attachments

        """)]
    for name, _files in sorted(files.items()):
        if not isinstance(_files, list):
            _files = [_files]
        links = []
        for file in _files:
            data = data_uri(file)
            links.append(':raw-html:`<a href="{data}" download="{filename}" draggable="true">{filename}</a>`'.format(
                data=data, filename=os.path.basename(file)))
        links = "\n\n          ".join(links)
        attachments.append('''
   .. container::
      :name: {name}

      {name}:
          {links}
            '''.format(name=name,
                       links=links))

    text = definitions + text + "\n\n" + "\n\n".join(attachments) + metadata

    overrides = dict()
    if template is not None:
        overrides["template"] = template
    if stylesheet is not None:
        overrides["stylesheet_path"] = stylesheet
    html = open(path, "w")
    publish_file(source=io.StringIO(text),
                 destination=html,
                 writer_name="html",
                 settings_overrides=overrides)
