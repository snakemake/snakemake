__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from asyncio import TaskGroup
from dataclasses import InitVar, dataclass, field
import os
import sys
import mimetypes
import base64
import textwrap
import datetime
import io
from typing import List, Optional
import uuid
import itertools
from collections import defaultdict
import hashlib
from pathlib import Path
import numbers


from docutils.parsers.rst.directives.images import Image, Figure
from docutils.parsers.rst import directives
from docutils.core import publish_file, publish_parts
from humanfriendly import format_size

from snakemake import script, wrapper, notebook
from snakemake.io.fmt import fmt_iofile
from snakemake.jobs import Job
from snakemake.report.common import data_uri_from_file, mime_from_file
from snakemake.rules import Rule
from snakemake.utils import format
from snakemake.logging import logger
from snakemake.io import (
    is_callable,
    is_flagged,
    get_flag_value,
    glob_wildcards,
    Wildcards,
    apply_wildcards,
    contains_wildcard,
)
from snakemake.exceptions import InputFunctionException, WorkflowError
from snakemake.script import Snakemake
from snakemake.common import (
    get_input_function_aux_params,
)
from snakemake import logging
from snakemake_interface_report_plugins.registry.plugin import Plugin as ReportPlugin
from snakemake_interface_report_plugins.settings import ReportSettingsBase
from snakemake_interface_report_plugins.interfaces import (
    CategoryInterface,
    RuleRecordInterface,
    ConfigFileRecordInterface,
    JobRecordInterface,
    FileRecordInterface,
)
from snakemake.common import get_report_id
from snakemake.exceptions import WorkflowError


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


# legacy report code
def report(
    text,
    path,
    stylesheet=None,
    defaultenc="utf8",
    template=None,
    metadata=None,
    **files,
):
    if stylesheet is None:
        os.path.join(os.path.dirname(__file__), "report.css")
    outmime, _ = mimetypes.guess_type(path)
    if outmime != "text/html":
        raise ValueError("Path to report output has to be an HTML file.")
    definitions = textwrap.dedent(
        """[]
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


async def expand_report_argument(item, wildcards, job):
    if is_callable(item):
        aux_params = get_input_function_aux_params(
            item, {"params": job.params, "input": job.input, "output": job.output}
        )
        io_items = ["input", "output"]
        if any(io_item in aux_params for io_item in io_items):
            # retrieve all input or output files from storage before evaluating function
            async with TaskGroup() as tg:
                for io_item in io_items:
                    if io_item in aux_params:
                        for f in aux_params[io_item]:
                            if f.is_storage:
                                tg.create_task(f.retrieve_from_storage())

        try:
            item = item(wildcards, **aux_params)
        except Exception as e:
            raise InputFunctionException(e, rule=job.rule, wildcards=wildcards)
    if isinstance(item, str):
        try:
            return apply_wildcards(item, wildcards)
        except AttributeError as e:
            raise WorkflowError("Failed to resolve wildcards.", e, rule=job.rule)
    elif isinstance(item, numbers.Number) and not isinstance(item, (int, float)):
        return str(item)
    else:
        return item


@dataclass(slots=True)
class Category(CategoryInterface):
    name: Optional[str]
    is_other: bool = field(init=False)
    id: str = field(init=False)

    def __post_init__(self):
        if self.name is None:
            self.name = "Other"

        self.is_other = self.name == "Other"
        h = hashlib.sha256()
        h.update(self.name.encode())
        self.id = h.hexdigest()

    def __eq__(self, other):
        return self.name.__eq__(other.name)

    def __hash__(self):
        return self.name.__hash__()

    def __lt__(self, other):
        if self.name == "other":
            return False
        return self.name.__lt__(other.name)


def render_iofile(iofile):
    if is_callable(iofile):
        return "<function>"
    else:
        return str(iofile)


@dataclass(slots=True)
class RuleRecord(RuleRecordInterface):
    job: InitVar
    job_rec: InitVar
    name: str = field(init=False)
    container_img_url: Optional[str] = field(init=False)
    conda_env: Optional[str] = field(init=False)
    n_jobs: int = field(init=False)
    id: str = field(init=False)
    language: str = field(init=False)
    source: str = field(init=False)

    def __post_init__(self, job, job_rec):
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

    def init_source(self):
        sources, language = None, None
        if self._rule.shellcmd is not None:
            sources = [self._rule.shellcmd]
            language = "bash"
        elif self._rule.script is not None and not contains_wildcard(self._rule.script):
            logger.info(f"Loading script code for rule {self.name}")
            _, source, language, _, _ = script.get_source(
                self._rule.script, self._rule.workflow.sourcecache, self._rule.basedir
            )
            sources = [source]
        elif self._rule.wrapper is not None and not contains_wildcard(
            self._rule.wrapper
        ):
            logger.info(f"Loading wrapper code for rule {self.name}")
            _, source, language, _, _ = script.get_source(
                wrapper.get_script(
                    self._rule.wrapper,
                    self._rule.workflow.sourcecache,
                    prefix=self._rule.workflow.workflow_settings.wrapper_prefix,
                ),
                self._rule.workflow.sourcecache,
            )
            sources = [source]
        elif self._rule.notebook is not None and not contains_wildcard(
            self._rule.notebook
        ):
            _, source, language, _, _ = script.get_source(
                self._rule.notebook, self._rule.workflow.sourcecache, self._rule.basedir
            )
            language = language.split("_")[1]
            sources = notebook.get_cell_sources(source)
        else:
            # A run directive. There is no easy way yet to obtain
            # the actual uncompiled source code.
            sources = []
            language = "python"

        self.source = "\n\n".join(sources)
        self.language = language

    def add(self, job_rec):
        self.n_jobs += 1

    @property
    def output(self):
        return [render_iofile(f) for f in self._rule.output]

    @property
    def input(self):
        return [render_iofile(f) for f in self._rule.input]

    def __eq__(self, other):
        return (
            self.name == other.name
            and self.conda_env == other.conda_env
            and self.container_img_url == other.container_img_url
        )


@dataclass(slots=True)
class ConfigfileRecord(ConfigFileRecordInterface):
    configfile: InitVar
    path: Path = field(init=False)
    source: str = field(init=False)

    def __post_init__(self, configfile):
        self.path = Path(configfile)
        self.source = self.path.read_text()


@dataclass(slots=True)
class JobRecord(JobRecordInterface):
    job: Job = field(init=False)
    rule: Rule = field(init=False)
    starttime: int = sys.maxsize
    endtime: int = 0
    output: list = field(default_factory=list)
    conda_env_file: Optional[Path] = field(init=False)
    container_img_url: Optional[Path] = field(init=False)


@dataclass(slots=True)
class FileRecord(FileRecordInterface):
    path: Path
    job: Job
    parent_path: Optional[Path] = None
    category: Optional[str] = None
    wildcards_overwrite: Optional[Wildcards] = None
    labels: Optional[dict] = None
    raw_caption: Optional[Path] = None
    aux_files: List[Path] = field(default_factory=list)
    name_overwrite: Optional[str] = None
    size: int = field(init=False)
    params: str = field(init=False)
    wildcards: str = field(init=False)
    mime: str = field(init=False)
    caption: str = field(init=False)
    id: str = field(init=False)
    target: str = field(init=False)

    def __post_init__(self):
        self.target = str(self.path.name)
        self.size = os.path.getsize(self.path)
        logger.info(f"Adding {self.name} ({format_size(self.size)}).")
        self.mime, _ = mime_from_file(self.path)

        self.id = get_report_id(self.parent_path or self.path)
        self.wildcards = logging.format_wildcards(self.raw_wildcards)
        self.params = (
            logging.format_dict(self.job.params)
            .replace("\n", r"\n")
            .replace('"', r"\"")
        )
        self.aux_files = self.aux_files or []

    @property
    def raw_wildcards(self):
        return (
            self.job.wildcards
            if self.wildcards_overwrite is None
            else self.wildcards_overwrite
        )

    def render(self, env, rst_links, categories, files):
        if self.raw_caption is not None:
            job = self.job
            snakemake = Snakemake(
                job.input,
                job.output,
                job.params,
                self.raw_wildcards,
                job.threads,
                job.resources,
                job.log,
                job.rule.module_globals["config"],
                job.rule.name,
                None,
            )

            try:
                caption = (
                    self.workflow.sourcecache.open(self.raw_caption).read() + rst_links
                )
                caption = env.from_string(caption).render(
                    snakemake=snakemake, categories=categories, files=files
                )
                self.caption = publish_parts(caption, writer_name="html")["body"]
            except Exception as e:
                raise WorkflowError(
                    "Error loading caption file {} of output marked for report.".format(
                        self.raw_caption.get_path_or_uri()
                    ),
                    e,
                )
        else:
            self.caption = ""

    @property
    def name(self):
        if self.name_overwrite:
            return self.name_overwrite
        return os.path.basename(self.path)

    @property
    def filename(self):
        if self.parent_path is None:
            return os.path.basename(self.path)
        else:
            return str(self.path.relative_to(self.parent_path))

    @property
    def workflow(self):
        return self.job.rule.workflow


async def expand_labels(labels, wildcards, job):
    if labels is None:
        return None
    labels = await expand_report_argument(labels, wildcards, job)

    if labels is None:
        return None

    if not isinstance(labels, dict) or not all(
        isinstance(col, str) or isinstance(col, numbers.Number)
        for col in labels.values()
    ):
        if isinstance(labels, dict):
            label_types = {name: type(value) for name, value in labels.items()}
        else:
            label_types = type(labels)
        raise WorkflowError(
            "Expected dict of strings or numbers as labels argument given to report flag. "
            f"Received: {label_types}",
            rule=job.rule,
        )
    return {
        name: await expand_report_argument(col, wildcards, job)
        for name, col in labels.items()
    }


async def auto_report(
    dag,
    report_plugin: ReportPlugin,
    report_settings: ReportSettingsBase,
):
    try:
        from jinja2 import Environment, PackageLoader, UndefinedError
    except ImportError as e:
        raise WorkflowError(
            "Python package jinja2 must be installed to create reports."
        )

    logger.info("Creating report...")

    persistence = dag.workflow.persistence
    results = defaultdict(lambda: defaultdict(list))
    records = defaultdict(JobRecord)
    recorded_files = set()

    env = Environment(
        trim_blocks=True,
        lstrip_blocks=True,
    )

    for job in dag.jobs:
        for f in job.output:
            meta = persistence.metadata(f)
            if not meta:
                logger.warning(
                    "Missing metadata for file {}. Maybe metadata "
                    "was deleted or it was created using an older "
                    "version of Snakemake. This is a non critical "
                    "warning.".format(f)
                )
                continue

            def get_time(rectime, metatime, sel_func):
                if metatime is None:
                    return rectime
                return sel_func(metatime, rectime)

            try:
                job_hash = meta["job_hash"]
                rule = meta["rule"]
                job_rec = records[(job_hash, rule)]
                job_rec.rule = rule
                job_rec.job = job
                job_rec.starttime = get_time(job_rec.starttime, meta["starttime"], min)
                job_rec.endtime = get_time(job_rec.endtime, meta["endtime"], max)
                job_rec.conda_env_file = None
                job_rec.conda_env = meta["conda_env"]
                job_rec.container_img_url = meta["container_img_url"]
                job_rec.output.append(f)
            except KeyError as e:
                logger.warning(
                    "Metadata for file {} was created with a too "
                    "old Snakemake version.".format(f)
                )

        for f in itertools.chain(job.output, job.input):
            if is_flagged(f, "report") and f not in recorded_files:
                if not await f.exists():
                    raise WorkflowError(
                        "File {} marked for report but does not exist.".format(f)
                    )
                report_obj = get_flag_value(f, "report")

                async def register_file(
                    f,
                    parent_path=None,
                    wildcards_overwrite=None,
                    aux_files=None,
                    name_overwrite=None,
                ):
                    wildcards = wildcards_overwrite or job.wildcards

                    async def expand_cat_name(cat_name, wildcards, job):
                        if cat_name is not None:
                            return await expand_report_argument(
                                cat_name, wildcards, job
                            )
                        else:
                            return cat_name

                    category = Category(
                        name=await expand_cat_name(report_obj.category, wildcards, job),
                    )
                    subcategory = Category(
                        name=await expand_cat_name(
                            report_obj.subcategory, wildcards, job
                        ),
                    )
                    labels = await expand_labels(report_obj.labels, wildcards, job)

                    results[category][subcategory].append(
                        FileRecord(
                            path=Path(f),
                            parent_path=(
                                Path(parent_path) if parent_path is not None else None
                            ),
                            job=job,
                            category=category,
                            raw_caption=report_obj.caption,
                            wildcards_overwrite=wildcards_overwrite,
                            aux_files=aux_files,
                            name_overwrite=name_overwrite,
                            labels=labels,
                        )
                    )
                    recorded_files.add(f)

                if f.is_storage:
                    await f.retrieve_from_storage()
                if os.path.isfile(f):
                    await register_file(f)
                elif os.path.isdir(f):
                    if report_obj.htmlindex:
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
                        await register_file(
                            os.path.join(f, report_obj.htmlindex),
                            aux_files=aux_files,
                            name_overwrite=f"{os.path.basename(f)}.html",
                        )
                    elif report_obj.patterns:
                        if not isinstance(report_obj.patterns, list):
                            raise WorkflowError(
                                "Invalid patterns given for report. Must be list.",
                                rule=job.rule,
                            )

                        found_something = False
                        for pattern in report_obj.patterns:
                            pattern = os.path.join(f, pattern)
                            wildcards = glob_wildcards(pattern)._asdict()
                            found_something |= len(wildcards) > 0
                            names = wildcards.keys()
                            for w in zip(*wildcards.values()):
                                w = dict(zip(names, w))
                                w.update(job.wildcards_dict)
                                w = Wildcards(fromdict=w)
                                subfile = apply_wildcards(pattern, w)
                                await register_file(
                                    subfile, parent_path=f, wildcards_overwrite=w
                                )
                        if not found_something:
                            logger.warning(
                                "No files found for patterns given to report marker "
                                f"in rule {job.rule.name} for output {fmt_iofile(f)}. Make sure "
                                "that the patterns are correctly specified."
                            )
                    else:
                        raise WorkflowError(
                            "Directory marked for report but neither file patterns "
                            "given via patterns=[...], nor htmlindex given. "
                            "See report documentation.",
                            rule=job.rule,
                        )

    for subcats in results.values():
        for catresults in subcats.values():
            catresults.sort(key=lambda res: res.name)

    # prepare per-rule information
    rules = defaultdict(list)
    for job_rec in records.values():
        rule = RuleRecord(job_rec.job, job_rec)
        if job_rec.rule not in rules:
            rules[job_rec.rule].append(rule)
            rule.init_source()
        else:
            merged = False
            for other in rules[job_rec.rule]:
                if rule == other:
                    other.add(job_rec)
                    merged = True
                    break
            if not merged:
                rules[job_rec.rule].append(rule)
                rule.init_source()
    # In theory there could be more than one rule with the same name kept from above.
    # For now, we just keep the first.
    rules = {rulename: items[0] for rulename, items in rules.items()}

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
    .. _{{ cat.name }}: javascript:app.showCategory('{{ cat.name|urlencode }}')
    {% endfor %}
    {% for res in files %}
    .. _{{ res.target }}: javascript:app.showResultInfo('{{ res.path|urlencode }}')
    {% endfor %}
    """
    )
    for cat, subcats in results.items():
        for subcat, catresults in subcats.items():
            for res in catresults:
                res.render(env, rst_links, results, files)

    # global description
    workflow_description = ""
    if dag.workflow.report_text:
        with dag.workflow.sourcecache.open(dag.workflow.report_text) as f:

            class Snakemake:
                config = dag.workflow.config

            workflow_description = f.read() + rst_links

            try:
                workflow_description = publish_parts(
                    env.from_string(workflow_description).render(
                        snakemake=Snakemake, categories=results, files=files
                    ),
                    writer_name="html",
                )["body"]
            except UndefinedError as e:
                raise WorkflowError(
                    "Error rendering global report caption {}:".format(
                        dag.workflow.report_text.get_path_or_uri()
                    ),
                    e,
                )

    reporter = report_plugin.reporter(
        rules,
        results,
        configfiles,
        sorted(records.values(), key=lambda rec: rec.rule),
        report_settings,
        workflow_description,
        dag=dag,
    )

    reporter.render()
    logger.info("Report created.")
