__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
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
from itertools import accumulate, chain
import urllib.parse
import hashlib
from zipfile import ZipFile
from pathlib import Path

import requests

from docutils.parsers.rst.directives.images import Image, Figure
from docutils.parsers.rst import directives
from docutils.core import publish_file, publish_parts

from snakemake import script, wrapper, notebook
from snakemake.utils import format
from snakemake.logging import logger
from snakemake.io import (
    is_flagged,
    get_flag_value,
    glob_wildcards,
    Wildcards,
    apply_wildcards,
    contains_wildcard,
)
from snakemake.exceptions import WorkflowError
from snakemake.script import Snakemake
from snakemake import __version__
from snakemake.common import num_if_possible, lazy_property
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
    stylesheet=None,
    defaultenc="utf8",
    template=None,
    metadata=None,
    **files
):
    if stylesheet is None:
        os.path.join(os.path.dirname(__file__), "report.css")
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
    def __init__(self, name, wildcards, job):
        if name is None:
            name = "Other"
            self.is_other = True
        else:
            self.is_other = False
            try:
                name = apply_wildcards(name, wildcards)
            except AttributeError as e:
                raise WorkflowError("Failed to resolve wildcards.", e, rule=job.rule)
        self.name = name
        h = hashlib.sha256()
        h.update(name.encode())
        self.id = h.hexdigest()

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
        self.container_img_url = job_rec.container_img_url
        self.conda_env = None
        self._conda_env_raw = None
        if job_rec.conda_env:
            self._conda_env_raw = base64.b64decode(job_rec.conda_env).decode()
            self.conda_env = yaml.load(self._conda_env_raw, Loader=yaml.Loader)
        self.n_jobs = 1
        self.id = uuid.uuid4()

    @lazy_property
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
        sources, language = None, None
        if self._rule.shellcmd is not None:
            sources = [self._rule.shellcmd]
            language = "bash"
        elif self._rule.script is not None and not contains_wildcard(self._rule.script):
            logger.info("Loading script code for rule {}".format(self.name))
            _, source, language = script.get_source(
                self._rule.script, self._rule.basedir
            )
            sources = [source.decode()]
        elif self._rule.wrapper is not None and not contains_wildcard(
            self._rule.wrapper
        ):
            logger.info("Loading wrapper code for rule {}".format(self.name))
            _, source, language = script.get_source(
                wrapper.get_script(
                    self._rule.wrapper, prefix=self._rule.workflow.wrapper_prefix
                )
            )
            sources = [source.decode()]
        elif self._rule.notebook is not None and not contains_wildcard(
            self._rule.notebook
        ):
            _, source, language = script.get_source(
                self._rule.notebook, self._rule.basedir
            )
            language = language.split("_")[1]
            sources = notebook.get_cell_sources(source)
        else:
            # A run directive. There is no easy way yet to obtain
            # the actual uncompiled source code.
            sources = []
            language = "python"

        try:
            lexer = get_lexer_by_name(language)

            highlighted = [
                highlight(
                    source,
                    lexer,
                    HtmlFormatter(linenos=True, cssclass="source", wrapcode=True),
                )
                for source in sources
            ]

            return highlighted
        except pygments.util.ClassNotFound:
            return [
                '<pre class="source"><code>{}</code></pre>'.format(source)
                for source in sources
            ]

    def add(self, job_rec):
        self.n_jobs += 1

    @property
    def output(self):
        return self._rule.output

    @property
    def input(self):
        return self._rule.input

    def __eq__(self, other):
        return (
            self.name == other.name
            and self.conda_env == other.conda_env
            and self.container_img_url == other.container_img_url
        )


class ConfigfileRecord:
    def __init__(self, configfile):
        self.path = Path(configfile)

    def code(self):
        try:
            from pygments.lexers import get_lexer_by_name
            from pygments.formatters import HtmlFormatter
            from pygments import highlight
        except ImportError:
            raise WorkflowError(
                "Python package pygments must be installed to create reports."
            )

        file_ext = self.path.suffix
        if file_ext in (".yml", ".yaml"):
            language = "yaml"
        elif file_ext == ".json":
            language = "json"
        else:
            raise ValueError(
                "Config file extension {} is not supported - must be YAML or JSON".format(
                    file_ext
                )
            )

        lexer = get_lexer_by_name(language)
        return highlight(
            self.path.read_text(),
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
        self.container_img_url = None


class FileRecord:
    def __init__(
        self,
        path,
        job,
        caption,
        env,
        category,
        wildcards_overwrite=None,
        mode_embedded=True,
        aux_files=None,
    ):
        self.mode_embedded = mode_embedded
        self.path = path
        self.target = os.path.basename(path)
        self.size = os.path.getsize(self.path)
        logger.info("Adding {} ({:.2g} MB).".format(self.name, self.size / 1e6))
        self.raw_caption = caption
        self.mime, _ = mime_from_file(self.path)

        h = hashlib.sha256()
        h.update(path.encode())

        self.id = h.hexdigest()
        self.job = job
        self._wildcards = (
            job.wildcards if wildcards_overwrite is None else wildcards_overwrite
        )
        self.wildcards = logging.format_wildcards(self._wildcards)
        self.params = (
            logging.format_dict(job.params).replace("\n", r"\n").replace('"', r"\"")
        )
        self.category = category

        self.aux_files = aux_files or []

        self.data_uri = self._data_uri()
        self.png_uri = self._png_uri()

    @lazy_property
    def png_content(self):
        assert self.is_img

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
                    ["convert", "-resize", max_spec, self.path, "png:-"], stderr=sp.PIPE
                )
                return png
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

    def _png_uri(self):
        if self.is_img:
            png = self.png_content
            if self.mode_embedded:
                if png is not None:
                    uri = data_uri(
                        png, os.path.basename(self.path) + ".png", mime="image/png"
                    )
                    return uri
            else:
                if png is not None:
                    return os.path.join("data/thumbnails", self.id)

    def _data_uri(self):
        if self.mode_embedded:
            return data_uri_from_file(self.path)
        else:
            return os.path.join("data/raw", self.id, self.name)

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
                self._wildcards,
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
                self.caption = json.dumps(
                    publish_parts(caption, writer_name="html")["body"]
                )
            except Exception as e:
                raise WorkflowError(
                    "Error loading caption file of output marked for report.", e
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


def auto_report(dag, path, stylesheet=None):
    try:
        from jinja2 import Template, Environment, PackageLoader
    except ImportError as e:
        raise WorkflowError(
            "Python package jinja2 must be installed to create reports."
        )

    mode_embedded = True
    if path.endswith(".zip"):
        mode_embedded = False
    elif not path.endswith(".html"):
        raise WorkflowError("Report file does not end with .html or .zip")

    custom_stylesheet = None
    if stylesheet is not None:
        try:
            with open(stylesheet) as s:
                custom_stylesheet = s.read()
        except (Exception, BaseException) as e:
            raise WorkflowError("Unable to read custom report stylesheet.", e)

    logger.info("Creating report...")

    env = Environment(
        loader=PackageLoader("snakemake", "report"),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["get_resource_as_string"] = get_resource_as_string

    persistence = dag.workflow.persistence
    results = defaultdict(lambda: defaultdict(list))
    records = defaultdict(JobRecord)
    recorded_files = set()
    for job in dag.jobs:
        for f in itertools.chain(job.expanded_output, job.input):
            if is_flagged(f, "report") and f not in recorded_files:
                if not f.exists:
                    raise WorkflowError(
                        "File {} marked for report but does " "not exist.".format(f)
                    )
                report_obj = get_flag_value(f, "report")

                def register_file(f, wildcards_overwrite=None, aux_files=None):
                    wildcards = wildcards_overwrite or job.wildcards
                    category = Category(
                        report_obj.category, wildcards=wildcards, job=job
                    )
                    subcategory = Category(
                        report_obj.subcategory, wildcards=wildcards, job=job
                    )

                    results[category][subcategory].append(
                        FileRecord(
                            f,
                            job,
                            report_obj.caption,
                            env,
                            category,
                            wildcards_overwrite=wildcards_overwrite,
                            mode_embedded=mode_embedded,
                            aux_files=aux_files,
                        )
                    )
                    recorded_files.add(f)

                if os.path.isfile(f):
                    register_file(f)
                elif os.path.isdir(f):
                    if report_obj.htmlindex:
                        if mode_embedded:
                            raise WorkflowError(
                                "Directory marked for report specifies htmlindex. "
                                "This is unsupported when requesting a pure HTML report. "
                                "Please use store as zip instead (--report report.zip)."
                            )
                        aux_files = []
                        index_found = False
                        for root, dirs, files in os.walk(f):
                            for name in files:
                                if name != ".snakemake_timestamp":
                                    filepath = os.path.join(root, name)
                                    if (
                                        os.path.relpath(filepath, f)
                                        != report_obj.htmlindex
                                    ):
                                        aux_files.append(filepath)
                                    else:
                                        index_found = True
                        if not index_found:
                            raise WorkflowError(
                                "Given htmlindex {} not found in directory "
                                "marked for report".format(report_obj.htmlindex)
                            )
                        register_file(
                            os.path.join(f, report_obj.htmlindex), aux_files=aux_files
                        )
                    elif report_obj.patterns:
                        if not isinstance(report_obj.patterns, list):
                            raise WorkflowError(
                                "Invalid patterns given for report. Must be list.",
                                rule=job.rule,
                            )

                        for pattern in report_obj.patterns:
                            pattern = os.path.join(f, pattern)
                            wildcards = glob_wildcards(pattern)._asdict()
                            names = wildcards.keys()
                            for w in zip(*wildcards.values()):
                                w = dict(zip(names, w))
                                w.update(job.wildcards_dict)
                                w = Wildcards(fromdict=w)
                                f = apply_wildcards(pattern, w)
                                register_file(f, wildcards_overwrite=w)
                    else:
                        raise WorkflowError(
                            "Directory marked for report but neither file patterns "
                            "given via patterns=[...], nor htmlindex given. "
                            "See report documentation.",
                            rule=job.rule,
                        )

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
                rec.container_img_url = meta["container_img_url"]
                rec.output.append(f)
            except KeyError as e:
                print(e)
                logger.warning(
                    "Metadata for file {} was created with a too "
                    "old Snakemake version.".format(f)
                )

    for subcats in results.values():
        for catresults in subcats.values():
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
        for subcat in cat.values()
        for res in subcat
        if res.target not in seen
    ]

    rst_links = textwrap.dedent(
        """

    .. _Workflow: javascript:show_panel('workflow')
    .. _Statistics: javascript:show_panel('statistics')
    {% for cat, catresults in categories|dictsort %}
    .. _{{ cat.name }}: javascript:show_panel("{{ cat.id }}")
    {% endfor %}
    {% for res in files %}
    .. _{{ res.target }}: javascript:show_panel("{{ res.category.id }}")
    {% endfor %}
    """
    )
    for cat, subcats in results.items():
        for subcat, catresults in subcats.items():
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
    results_size = sum(
        res.size
        for cat in results.values()
        for subcat in cat.values()
        for res in subcat
    )

    try:
        from pygments.formatters import HtmlFormatter
    except ImportError:
        raise WorkflowError(
            "Python package pygments must be installed to create reports."
        )

    template = env.get_template("report.html.jinja2")

    logger.info("Downloading resources and rendering HTML.")

    rendered = template.render(
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
        custom_stylesheet=custom_stylesheet,
        mode_embedded=mode_embedded,
    )

    # TODO look into supporting .WARC format, also see (https://webrecorder.io)

    if not mode_embedded:
        with ZipFile(path, mode="w") as zipout:
            folder = Path(Path(path).stem)
            # store results in data folder
            for subcats in results.values():
                for catresults in subcats.values():
                    for result in catresults:
                        # write raw data
                        zipout.write(result.path, str(folder.joinpath(result.data_uri)))
                        # write thumbnail
                        if result.is_img and result.png_content:
                            zipout.writestr(
                                str(folder.joinpath(result.png_uri)), result.png_content
                            )
                        # write aux files
                        parent = folder.joinpath(result.data_uri).parent
                        for aux_path in result.aux_files:
                            # print(aux_path, parent, str(parent.joinpath(os.path.relpath(aux_path, os.path.dirname(result.path)))))
                            zipout.write(
                                aux_path,
                                str(
                                    parent.joinpath(
                                        os.path.relpath(
                                            aux_path, os.path.dirname(result.path)
                                        )
                                    )
                                ),
                            )

            # write report html
            zipout.writestr(str(folder.joinpath("report.html")), rendered)
    else:
        with open(path, "w", encoding="utf-8") as htmlout:
            htmlout.write(rendered)

    logger.info("Report created: {}.".format(path))
