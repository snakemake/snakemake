__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import sys
import mimetypes
import base64
import textwrap
import datetime
import io
import uuid
import json
import time
import shutil
import subprocess as sp
import itertools
import csv
from collections import namedtuple, defaultdict

import requests

from docutils.parsers.rst.directives.images import Image, Figure
from docutils.parsers.rst import directives
from docutils.core import publish_file, publish_parts

from snakemake import script, wrapper
from snakemake.utils import format
from snakemake.logging import logger
from snakemake.io import is_flagged, get_flag_value
from snakemake.exceptions import WorkflowError
from snakemake.script import Snakemake
from snakemake import __version__
from snakemake.common import num_if_possible
from snakemake import logging


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
        self.options["uri"] = data_uri_from_file(reference)
        return result


# Create (and register) new image:: and figure:: directives that use a base64
# data URI instead of pointing to a filename.


class EmbeddedImage(Image, EmbeddedMixin):
    pass


directives.register_directive("embeddedimage", EmbeddedImage)


class EmbeddedFigure(Figure, EmbeddedMixin):
    pass


directives.register_directive("embeddedfigure", EmbeddedFigure)


def data_uri(data, filename, encoding="utf8", mime="text/plain"):
    """Craft a base64 data URI from file with proper encoding and mimetype."""
    data = base64.b64encode(data)
    uri = "data:{mime};charset={charset};filename={filename};base64,{data}" "".format(
        filename=filename, mime=mime, charset=encoding, data=data.decode("utf-8")
    )
    return uri


def mime_from_file(file):
    mime, encoding = mimetypes.guess_type(file)
    if mime is None:
        mime = "text/plain"
        logger.info(
            "Could not detect mimetype for {}, assuming " "text/plain.".format(file)
        )
    return mime, encoding


def data_uri_from_file(file, defaultenc="utf8"):
    """Craft a base64 data URI from file with proper encoding and mimetype."""
    mime, encoding = mime_from_file(file)
    if encoding is None:
        encoding = defaultenc
    with open(file, "rb") as f:
        return data_uri(f.read(), os.path.basename(file), encoding, mime)


def report(
    text,
    path,
    stylesheet=os.path.join(os.path.dirname(__file__), "report.css"),
    defaultenc="utf8",
    template=None,
    metadata=None,
    **files
):
    outmime, _ = mimetypes.guess_type(path)
    if outmime != "text/html":
        raise ValueError("Path to report output has to be an HTML file.")
    definitions = textwrap.dedent(
        """
    .. role:: raw-html(raw)
       :format: html

    """
    )

    metadata = textwrap.dedent(
        """

    .. container::
       :name: metadata

       {metadata}{date}

    """
    ).format(
        metadata=metadata + " | " if metadata else "",
        date=datetime.date.today().isoformat(),
    )

    text = format(textwrap.dedent(text), stepout=3)

    attachments = []
    if files:
        attachments = [
            textwrap.dedent(
                """
            .. container::
               :name: attachments

            """
            )
        ]
        for name, _files in sorted(files.items()):
            if not isinstance(_files, list):
                _files = [_files]
            links = []
            for file in sorted(_files):
                data = data_uri_from_file(file)
                links.append(
                    ':raw-html:`<a href="{data}" download="{filename}" draggable="true">{filename}</a>`'.format(
                        data=data, filename=os.path.basename(file)
                    )
                )
            links = "\n\n              ".join(links)
            attachments.append(
                """
       .. container::
          :name: {name}

          {name}:
              {links}
                """.format(
                    name=name, links=links
                )
            )

    text = definitions + text + "\n\n" + "\n\n".join(attachments) + metadata

    overrides = dict()
    if template is not None:
        overrides["template"] = template
    if stylesheet is not None:
        overrides["stylesheet_path"] = stylesheet
    html = open(path, "w")
    publish_file(
        source=io.StringIO(text),
        destination=html,
        writer_name="html",
        settings_overrides=overrides,
    )


class Category:
    def __init__(self, name):
        if name is None:
            name = "Other"
        self.name = name
        self.id = "results-{name}".format(name=name.replace(" ", "_"))
        self.content_id = self.id + "-content"

    def __eq__(self, other):
        return self.name.__eq__(other.name)

    def __hash__(self):
        return self.name.__hash__()

    def __lt__(self, other):
        if self.name == "other":
            return False
        return self.name.__lt__(other.name)


class RuleRecord:
    def __init__(self, job, job_rec):
        import yaml

        self.name = job_rec.rule
        self._rule = job.rule
        self.singularity_img_url = job_rec.singularity_img_url
        self.conda_env = None
        self._conda_env_raw = None
        if job_rec.conda_env:
            self._conda_env_raw = base64.b64decode(job_rec.conda_env).decode()
            self.conda_env = yaml.load(self._conda_env_raw, Loader=yaml.Loader)
        self.n_jobs = 1
        self.output = list(job_rec.output)
        self.id = uuid.uuid4()

    def code(self):
        try:
            from pygments.lexers import get_lexer_by_name
            from pygments.formatters import HtmlFormatter
            from pygments import highlight
            import pygments.util
        except ImportError:
            raise WorkflowError(
                "Python package pygments must be installed to create reports."
            )
        source, language = None, None
        if self._rule.shellcmd is not None:
            source = self._rule.shellcmd
            language = "bash"
        elif self._rule.script is not None:
            logger.info("Loading script code for rule {}".format(self.name))
            _, source, language = script.get_source(
                self._rule.script, self._rule.basedir
            )
            source = source.decode()
        elif self._rule.wrapper is not None:
            logger.info("Loading wrapper code for rule {}".format(self.name))
            _, source, language = script.get_source(
                wrapper.get_script(
                    self._rule.wrapper, prefix=self._rule.workflow.wrapper_prefix
                )
            )
            source = source.decode()

        try:
            lexer = get_lexer_by_name(language)
            return highlight(
                source,
                lexer,
                HtmlFormatter(linenos=True, cssclass="source", wrapcode=True),
            )
        except pygments.util.ClassNotFound:
            return "<pre><code>source</code></pre>"

    def add(self, job_rec):
        self.n_jobs += 1
        self.output.extend(job_rec.output)

    def __eq__(self, other):
        return (
            self.name == other.name
            and self.conda_env == other.conda_env
            and self.singularity_img_url == other.singularity_img_url
        )


class ConfigfileRecord:
    def __init__(self, configfile):
        self.name = configfile

    def code(self):
        try:
            from pygments.lexers import get_lexer_by_name
            from pygments.formatters import HtmlFormatter
            from pygments import highlight
        except ImportError:
            raise WorkflowError(
                "Python package pygments must be installed to create reports."
            )

        language = (
            "yaml"
            if self.name.endswith(".yaml") or self.name.endswith(".yml")
            else "json"
        )
        lexer = get_lexer_by_name(language)
        with open(self.name) as f:
            return highlight(
                f.read(),
                lexer,
                HtmlFormatter(linenos=True, cssclass="source", wrapcode=True),
            )


class JobRecord:
    def __init__(self):
        self.job = None
        self.rule = None
        self.starttime = sys.maxsize
        self.endtime = 0
        self.output = []
        self.conda_env_file = None
        self.singularity_img_url = None


class FileRecord:
    def __init__(self, path, job, caption, env, category):
        self.path = path
        self.target = os.path.basename(path)
        self.size = os.path.getsize(self.path)
        logger.info("Adding {} ({:.2g} MB).".format(self.name, self.size / 1e6))
        self.raw_caption = caption
        self.mime, _ = mime_from_file(self.path)
        self.id = uuid.uuid4()
        self.job = job
        self.wildcards = logging.format_wildcards(job.wildcards)
        self.params = logging.format_dict(job.params)
        self.png_uri = None
        self.category = category
        if self.is_img:
            convert = shutil.which("convert")
            if convert is not None:
                try:
                    # 2048 aims at a reasonable balance between what displays
                    # can show in a png-preview image and what renders quick
                    # into a small enough png
                    max_width = "2048"
                    max_height = "2048"
                    # '>' means only larger images scaled down to within max-dimensions
                    max_spec = max_width + "x" + max_height + ">"
                    png = sp.check_output(
                        ["convert", "-resize", max_spec, self.path, "png:-"],
                        stderr=sp.PIPE,
                    )
                    uri = data_uri(
                        png, os.path.basename(self.path) + ".png", mime="image/png"
                    )
                    self.png_uri = uri
                except sp.CalledProcessError as e:
                    logger.warning(
                        "Failed to convert image to png with "
                        "imagemagick convert: {}".format(e.stderr)
                    )
            else:
                logger.warning(
                    "Command convert not in $PATH. Install "
                    "imagemagick in order to have embedded "
                    "images and pdfs in the report."
                )
        if self.is_table:
            if self.size > 1e6:
                logger.warning(
                    "Table {} >1MB. Rendering as generic file.".format(self.path)
                )
            else:
                with open(self.path) as table:
                    dialect = None
                    for prefix in range(10, 17):
                        try:
                            table.seek(0)
                            dialect = csv.Sniffer().sniff(table.read(prefix))
                            break
                        except csv.Error:
                            pass
                        except UnicodeDecodeError:
                            # table is not readable as UTF-8
                            break
                    if dialect is None:
                        logger.warning(
                            "Failed to infer CSV/TSV dialect from table {}. "
                            "Rendering as generic file.".format(self.path)
                        )
                    else:
                        table.seek(0)
                        reader = csv.reader(table, dialect)
                        columns = next(reader)
                        table = map(lambda row: list(map(num_if_possible, row)), reader)
                        template = env.get_template("table.html")
                        html = template.render(
                            columns=columns, table=table, name=self.name
                        ).encode()
                        self.mime = "text/html"
                        self.path = os.path.basename(self.path) + ".html"
                        self.data_uri = data_uri(html, self.path, mime=self.mime)
                        return
        # fallback
        self.data_uri = data_uri_from_file(path)

    def render(self, env, rst_links, categories, files):
        if self.raw_caption is not None:
            try:
                from jinja2 import Template
            except ImportError as e:
                raise WorkflowError(
                    "Python package jinja2 must be installed to create reports."
                )

            job = self.job
            snakemake = Snakemake(
                job.input,
                job.output,
                job.params,
                job.wildcards,
                job.threads,
                job.resources,
                job.log,
                job.dag.workflow.config,
                job.rule.name,
                None,
            )

            try:
                caption = open(self.raw_caption).read() + rst_links
                caption = env.from_string(caption).render(
                    snakemake=snakemake, categories=categories, files=files
                )
                self.caption = publish_parts(caption, writer_name="html")["body"]
            except Exception as e:
                raise WorkflowError(
                    "Error loading caption file of output " "marked for report.", e
                )

    @property
    def is_img(self):
        web_safe = {
            "image/gif",
            "image/jpeg",
            "image/png",
            "image/svg+xml",
            "application/pdf",
        }
        return self.mime in web_safe

    @property
    def is_text(self):
        return self.is_table or self.mime == "text/plain"

    @property
    def is_table(self):
        return self.mime in {"text/csv", "text/tab-separated-values"}

    @property
    def is_vega(self):
        return (
            self.mime == "application/json"
            and self.path.endswith(".vl.json")
            or self.path.endswith(".vg.json")
        )

    @property
    def icon(self):
        if self.is_img:
            return "image"
        elif self.is_text:
            return "file-text"
        else:
            return "file"

    @property
    def name(self):
        return os.path.basename(self.path)


def rulegraph_d3_spec(dag):
    try:
        import networkx as nx
        from networkx.drawing.nx_agraph import graphviz_layout
        from networkx.readwrite import json_graph
    except ImportError as e:
        raise WorkflowError(
            "Python packages networkx and pygraphviz must be "
            "installed to create reports.",
            e,
        )

    g = nx.DiGraph()
    g.add_nodes_from(sorted(job.rule.name for job in dag.jobs))

    for job in dag.jobs:
        target = job.rule.name
        for dep in dag.dependencies[job]:
            source = dep.rule.name
            g.add_edge(source, target)

    pos = graphviz_layout(g, "dot", args="-Grankdir=BT")
    xmax = max(x for x, y in pos.values()) + 100  # add offset to account for labels
    ymax = max(y for x, y in pos.values())

    def encode_node(node):
        x, y = pos[node]
        return {"rule": node, "fx": x, "fy": y}

    nodes = list(map(encode_node, g.nodes))
    idx = {node: i for i, node in enumerate(g.nodes)}
    links = [{"target": idx[u], "source": idx[v], "value": 1} for u, v in g.edges]
    return {"nodes": nodes, "links": links}, xmax, ymax


def get_resource_as_string(url):
    r = requests.get(url)
    if r.status_code == requests.codes.ok:
        return r.text
    raise WorkflowError(
        "Failed to download resource needed for " "report: {}".format(url)
    )


def auto_report(dag, path):
    try:
        from jinja2 import Template, Environment, PackageLoader
    except ImportError as e:
        raise WorkflowError(
            "Python package jinja2 must be installed to create reports."
        )

    if not path.endswith(".html"):
        raise WorkflowError("Report file does not end with .html")

    logger.info("Creating report...")

    env = Environment(
        loader=PackageLoader("snakemake", "report"),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["get_resource_as_string"] = get_resource_as_string

    persistence = dag.workflow.persistence
    results = defaultdict(list)
    records = defaultdict(JobRecord)
    recorded_files = set()
    for job in dag.jobs:
        for f in itertools.chain(job.expanded_output, job.input):
            if is_flagged(f, "report") and f not in recorded_files:
                if not f.exists:
                    raise WorkflowError(
                        "File {} marked for report but does " "not exist.".format(f)
                    )
                if os.path.isfile(f):
                    report_obj = get_flag_value(f, "report")
                    category = Category(report_obj.category)
                    results[category].append(
                        FileRecord(f, job, report_obj.caption, env, category)
                    )
                    recorded_files.add(f)

        for f in job.expanded_output:
            meta = persistence.metadata(f)
            if not meta:
                logger.warning(
                    "Missing metadata for file {}. Maybe metadata "
                    "was deleted or it was created using an older "
                    "version of Snakemake. This is a non critical "
                    "warning.".format(f)
                )
                continue
            try:
                job_hash = meta["job_hash"]
                rule = meta["rule"]
                rec = records[(job_hash, rule)]
                rec.rule = rule
                rec.job = job
                rec.starttime = min(rec.starttime, meta["starttime"])
                rec.endtime = max(rec.endtime, meta["endtime"])
                rec.conda_env_file = None
                rec.conda_env = meta["conda_env"]
                rec.singularity_img_url = meta["singularity_img_url"]
                rec.output.append(f)
            except KeyError as e:
                print(e)
                logger.warning(
                    "Metadata for file {} was created with a too "
                    "old Snakemake version.".format(f)
                )

    for catresults in results.values():
        catresults.sort(key=lambda res: res.name)

    # prepare runtimes
    runtimes = [
        {"rule": rec.rule, "runtime": rec.endtime - rec.starttime}
        for rec in sorted(records.values(), key=lambda rec: rec.rule)
    ]

    # prepare end times
    timeline = [
        {
            "rule": rec.rule,
            "starttime": datetime.datetime.fromtimestamp(rec.starttime).isoformat(),
            "endtime": datetime.datetime.fromtimestamp(rec.endtime).isoformat(),
        }
        for rec in sorted(records.values(), key=lambda rec: rec.rule)
    ]

    # prepare per-rule information
    rules = defaultdict(list)
    for rec in records.values():
        rule = RuleRecord(rec.job, rec)
        if rec.rule not in rules:
            rules[rec.rule].append(rule)
        else:
            merged = False
            for other in rules[rec.rule]:
                if rule == other:
                    other.add(rec)
                    merged = True
                    break
            if not merged:
                rules[rec.rule].append(rule)

    # rulegraph
    rulegraph, xmax, ymax = rulegraph_d3_spec(dag)

    # configfiles
    configfiles = [ConfigfileRecord(f) for f in dag.workflow.configfiles]

    seen = set()
    files = [
        seen.add(res.target) or res
        for cat in results.values()
        for res in cat
        if res.target not in seen
    ]

    rst_links = textwrap.dedent(
        """

    .. _Results: #results
    .. _Rules: #rules
    .. _Statistics: #stats
    {% for cat, catresults in categories|dictsort %}
    .. _{{ cat.name }}: #{{ cat.id }}
    {% for res in files %}
    .. _{{ res.target }}: #{{ res.id }}
    {% endfor %}
    {% endfor %}
    .. _
    """
    )
    for cat, catresults in results.items():
        for res in catresults:
            res.render(env, rst_links, results, files)

    # global description
    text = ""
    if dag.workflow.report_text:
        with open(dag.workflow.report_text) as f:

            class Snakemake:
                config = dag.workflow.config

            text = f.read() + rst_links
            text = publish_parts(
                env.from_string(text).render(
                    snakemake=Snakemake, categories=results, files=files
                ),
                writer_name="html",
            )["body"]

    # record time
    now = "{} {}".format(datetime.datetime.now().ctime(), time.tzname[0])
    results_size = sum(res.size for cat in results.values() for res in cat)

    try:
        from pygments.formatters import HtmlFormatter
    except ImportError:
        raise WorkflowError(
            "Python package pygments must be installed to create reports."
        )

    # render HTML
    template = env.get_template("report.html")
    with open(path, "w", encoding="utf-8") as out:
        out.write(
            template.render(
                results=results,
                results_size=results_size,
                configfiles=configfiles,
                text=text,
                rulegraph_nodes=rulegraph["nodes"],
                rulegraph_links=rulegraph["links"],
                rulegraph_width=xmax + 20,
                rulegraph_height=ymax + 20,
                runtimes=runtimes,
                timeline=timeline,
                rules=[rec for recs in rules.values() for rec in recs],
                version=__version__,
                now=now,
                pygments_css=HtmlFormatter(style="trac").get_style_defs(".source"),
            )
        )
    logger.info("Report created.")
