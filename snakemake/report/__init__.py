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
from collections import namedtuple, defaultdict

import requests

from docutils.parsers.rst.directives.images import Image, Figure
from docutils.parsers.rst import directives
from docutils.core import publish_file, publish_parts

from snakemake.utils import format
from snakemake.logging import logger
from snakemake.io import is_flagged, get_flag_value
from snakemake.exceptions import WorkflowError
from snakemake.script import Snakemake
from snakemake import __version__


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
        self.options['uri'] = data_uri(reference)[0]
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
    return uri, mime


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

    attachments = []
    if files:
        attachments = [textwrap.dedent("""
            .. container::
               :name: attachments

            """)]
        for name, _files in sorted(files.items()):
            if not isinstance(_files, list):
                _files = [_files]
            links = []
            for file in sorted(_files):
                data, _ = data_uri(file)
                links.append(':raw-html:`<a href="{data}" download="{filename}" draggable="true">{filename}</a>`'.format(
                    data=data, filename=os.path.basename(file)))
            links = "\n\n              ".join(links)
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


class RuleRecord:
    def __init__(self, job, job_rec):
        import yaml
        self.name = job_rec.rule
        self.singularity_img_url = job_rec.singularity_img_url
        self.conda_env = None
        self._conda_env_raw = None
        if job_rec.conda_env:
            self._conda_env_raw = base64.b64decode(job_rec.conda_env).decode()
            self.conda_env = yaml.load(self._conda_env_raw)
        self.n_jobs = 1
        self.output = list(job_rec.output)
        self.id = uuid.uuid4()

    def add(self, job_rec):
        self.n_jobs += 1
        self.output.extend(job_rec.output)

    def __eq__(self, other):
        return (self.name == other.name and
                self.conda_env == other.conda_env and
                self.singularity_img_url == other.singularity_img_url)


class JobRecord:
    def __init__(self):
        self.rule = None
        self.starttime = sys.maxsize
        self.endtime = 0
        self.output = []
        self.conda_env_file = None
        self.singularity_img_url = None



class FileRecord:
    def __init__(self, path, job, caption=None):
        self.caption = ""
        if caption is not None:
            try:
                from jinja2 import Template
            except ImportError as e:
                raise WorkflowError("Pyhton package jinja2 must be installed to create reports.")

            snakemake = Snakemake(job.input, job.output, job.params, job.wildcards,
                                  job.threads, job.resources, job.log,
                                  job.dag.workflow.config, job.rule.name)
            try:
                caption = Template(open(caption).read()).render(snakemake=snakemake)
                self.caption = publish_parts(caption, writer_name="html")["body"]
            except Exception as e:
                raise WorkflowError("Error loading caption file of output "
                                    "marked for report.", e)
        self.data, self.mime = data_uri(path)
        self.id = uuid.uuid4()
        self.path = path

    @property
    def is_img(self):
        web_safe = {"image/gif", "image/jpeg", "image/png", "image/svg+xml"}
        return self.mime in web_safe

    @property
    def name(self):
        return os.path.basename(self.path)

    @property
    def size(self):
        """Return size in MB."""
        return os.path.getsize(self.path) / 1024


def rulegraph_d3_spec(dag):
    try:
        import networkx as nx
        from networkx.drawing.nx_agraph import graphviz_layout
        from networkx.readwrite import json_graph
    except ImportError as e:
        raise WorkflowError("Python packages networkx and pygraphviz must be "
                            "installed to create reports.", e)

    g = nx.DiGraph()
    g.add_nodes_from(sorted(job.rule.name for job in dag.jobs))

    for job in dag.jobs:
        target = job.rule.name
        for dep in dag.dependencies[job]:
            source = dep.rule.name
            g.add_edge(source, target)

    pos = graphviz_layout(g, "dot", args="-Grankdir=LR")
    xmax = max(x for x, y in pos.values())
    ymax = max(y for x, y in pos.values())

    def encode_node(node):
        x, y = pos[node]
        return {"rule": node, "fx": x, "fy": y}

    nodes = list(map(encode_node, g.nodes))
    idx = {node: i for i, node in enumerate(g.nodes)}
    links = [{"target": idx[u], "source": idx[v], "value": 1}
             for u, v in g.edges]
    return {"nodes": nodes, "links": links}, xmax, ymax


def get_resource_as_string(url):
    r = requests.get(url)
    if r.status_code == requests.codes.ok:
        return r.text
    raise WorkflowError("Failed to download resource needed for "
                        "report: {}".format(url))


def auto_report(dag, path):
    try:
        from jinja2 import Template, Environment, PackageLoader
    except ImportError as e:
        raise WorkflowError("Pyhton package jinja2 must be installed to create reports.")

    if not path.endswith(".html"):
        raise WorkflowError("Report file does not end with .html")

    persistence = dag.workflow.persistence
    results = []
    records = defaultdict(JobRecord)
    for job in dag.jobs:
        for f in job.expanded_output:
            if is_flagged(f, "report"):
                if not f.exists:
                    raise WorkflowError("Output file {} marked for report but does "
                                        "not exist.")
                if os.path.isfile(f):
                    results.append(FileRecord(f, job, get_flag_value(f, "report")))

            meta = persistence.metadata(f)
            if not meta:
                logger.warning("Missing metadata for file {}".format(f))
                continue
            try:
                job_hash = meta["job_hash"]
                rule = meta["rule"]
                rec = records[(job_hash, rule)]
                rec.rule = rule
                rec.starttime = min(rec.starttime, meta["starttime"])
                rec.endtime = max(rec.endtime, meta["endtime"])
                rec.conda_env_file = None
                rec.conda_env = meta["conda_env"]
                rec.singularity_img_url = meta["singularity_img_url"]
                rec.output.append(f)
            except KeyError as e:
                print(e)
                logger.warning("Metadata for file {} was created with a too "
                               "old Snakemake version.".format(f))

    results.sort(key=lambda res: res.name)

    # prepare runtimes
    runtimes = [{"rule": rec.rule, "runtime": rec.endtime - rec.starttime}
                for rec in sorted(records.values(), key=lambda rec: rec.rule)]

    # prepare end times
    endtimes = [{"rule": rec.rule,
                 "endtime": datetime.datetime.fromtimestamp(rec.endtime).isoformat()}
                for rec in sorted(records.values(), key=lambda rec: rec.rule)]

    # prepare per-rule information
    rules = defaultdict(list)
    for rec in records.values():
        rule = RuleRecord(job, rec)
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

    # global description
    text = ""
    if dag.workflow.report_text:
        with open(dag.workflow.report_text) as f:
            class Snakemake:
                config = dag.workflow.config
            text = publish_parts(Template(f.read()).render(snakemake=Snakemake),
                                 writer_name="html")["body"]

    # rulegraph
    rulegraph, xmax, ymax = rulegraph_d3_spec(dag)

    # compose html
    env = Environment(loader=PackageLoader("snakemake", "report"))
    env.filters["get_resource_as_string"] = get_resource_as_string

    # record time
    now = "{} {}".format(datetime.datetime.now().ctime(), time.tzname[0])

    template = env.get_template("report.html")
    with open(path, "w") as out:
        out.write(template.render(results=results,
                                  results_size=sum(res.size for res in results),
                                  text=text,
                                  rulegraph_nodes=rulegraph["nodes"],
                                  rulegraph_links=rulegraph["links"],
                                  rulegraph_width=xmax + 20,
                                  rulegraph_height=ymax + 20,
                                  runtimes=runtimes,
                                  endtimes=endtimes,
                                  rules=[rec for recs in rules.values() for rec in recs],
                                  version=__version__.split("+")[0],
                                  now=now))
    logger.info("Report created.")
